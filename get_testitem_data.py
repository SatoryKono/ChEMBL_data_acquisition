"""Command line interface for retrieving ChEMBL and PubChem compound data."""

from __future__ import annotations

import argparse
import logging
from typing import Sequence

import pandas as pd
import requests
from library.config import Config, ensure_dirs, print_config

from library import chembl_library as cl
from library import pubchem_library as pl
from library import io
from library.cli import (
    apply_config_overrides,
    build_parser as base_parser,
    configure_logging,
)
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from pandera.errors import SchemaErrors
from schemas import TestitemsSchema, normalize_testitems

logger = logging.getLogger(__name__)


def add_pubchem_data(df: pd.DataFrame, cfg: pl.PubChemCfg) -> pd.DataFrame:
    """Augment ChEMBL records with PubChem information.

    For each canonical SMILES string in ``df``, the function looks up the
    corresponding PubChem CID and basic chemical properties. The PubChem
    fields are appended to the input frame. If a SMILES string cannot be
    resolved, empty values are inserted.

    Parameters
    ----------
    df:
        Data frame returned by :func:`library.chembl_library.get_testitem`.
    cfg:
        PubChem configuration options.

    Returns
    -------
    pandas.DataFrame
        ``df`` with additional PubChem columns.

    """
    if df.empty or "molecule_structures.canonical_smiles" not in df.columns:
        return df

    smiles_list = df["molecule_structures.canonical_smiles"].fillna("").tolist()
    # ``dict.fromkeys`` preserves the order of first occurrence while
    # removing duplicates. This allows progress output to reflect the
    # deterministic iteration order of SMILES strings.
    unique_smiles = [s for s in dict.fromkeys(smiles_list) if s]

    total = len(unique_smiles)
    if total:
        logger.info("Fetching PubChem data for %d unique SMILES", total)
    else:
        logger.info("No SMILES strings available for PubChem lookup")

    records: dict[str, dict[str, str]] = {}
    for idx, smi in enumerate(unique_smiles, start=1):
        logger.info("PubChem lookup %d/%d", idx, total)
        cid = pl.get_cid_from_smiles(smi, cfg) or ""
        first_cid = cid.split("|")[0] if cid else ""
        if first_cid:
            props = pl.get_properties(first_cid, cfg)
            records[smi] = {
                "pubchem_cid": first_cid,
                "pubchem_iupac_name": props.IUPACName,
                "pubchem_molecular_formula": props.MolecularFormula,
                "pubchem_isomeric_smiles": props.iSMILES,
                "pubchem_canonical_smiles": props.cSMILES,
                "pubchem_inchi": props.InChI,
                "pubchem_inchikey": props.InChIKey,
            }
        else:
            records[smi] = {
                "pubchem_cid": "",
                "pubchem_iupac_name": "",
                "pubchem_molecular_formula": "",
                "pubchem_isomeric_smiles": "",
                "pubchem_canonical_smiles": "",
                "pubchem_inchi": "",
                "pubchem_inchikey": "",
            }

    empty = {
        "pubchem_cid": "",
        "pubchem_iupac_name": "",
        "pubchem_molecular_formula": "",
        "pubchem_isomeric_smiles": "",
        "pubchem_canonical_smiles": "",
        "pubchem_inchi": "",
        "pubchem_inchikey": "",
    }
    pubchem_rows = [records.get(smi, empty) for smi in smiles_list]
    pubchem_df = pd.DataFrame(pubchem_rows)
    return pd.concat([df.reset_index(drop=True), pubchem_df], axis=1)


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute compound retrieval from the ChEMBL API and augment with PubChem data.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    try:
        ids = io.read_ids(
            args.input_csv,
            column=args.column,
            cfg=cfg.io,
            sep=args.sep,
            encoding=args.encoding,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    logger.info("Retrieved %d identifiers", len(ids))
    logger.info("Fetching ChEMBL data in chunks of %d", args.chunk_size)
    try:
        df = cl.get_testitem(
            ids, cfg=cfg.api, chunk_size=args.chunk_size, timeout=args.timeout
        )
    except (requests.RequestException, ValueError) as exc:
        logger.error("failed to retrieve compounds: %s", exc)
        return 1
    logger.info("Retrieved %d rows from ChEMBL", len(df))
    logger.info("Augmenting results with PubChem data")
    df = add_pubchem_data(df, cfg.pubchem)
    logger.info("PubChem augmentation completed")
    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    df = normalize_testitems(df)
    exit_code = 0
    try:
        df = TestitemsSchema.validate(df, lazy=True)
    except SchemaErrors as exc:
        failure_path = output.with_name(f"{output.stem}_failure_cases.csv")
        errors = SidecarErrors()
        for row in exc.failure_cases.to_dict("records"):
            errors.add_error(row)
        errors.save(failure_path)
        logger.error(
            "validation failed; wrote %d failure cases to %s",
            len(exc.failure_cases),
            failure_path,
        )
        df = exc.validated_data  # type: ignore[attr-defined]
        exit_code = 1
    try:
        io.write_csv(df, output, cfg=cfg, sep=args.sep, encoding=args.encoding)
        logger.info("Wrote %d rows to %s", len(df), output)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1
    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return exit_code


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""
    parser = base_parser(
        "ChEMBL and PubChem compound data utilities",
        column="molecule_chembl_id",
        chunk_size=5,
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    parser.set_defaults(func=run_chembl)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        cfg: Config = apply_config_overrides(
            args, parser, args.config, mapping={"timeout": "api.timeout_read"}
        )
        if args.print_config:
            print_config(cfg)
            return 0
        ensure_dirs(cfg)
        configure_logging(args.log_level, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error("%s", exc)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("failed to set up directories: %s", exc)
        return 1
    return args.func(cfg, args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
