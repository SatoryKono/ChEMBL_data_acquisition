"""Command line interface for retrieving ChEMBL and PubChem compound data."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

import pandas as pd

from library import chembl_library as cl
from library import pubchem_library as pl
from library import io

logger = logging.getLogger(__name__)


def add_pubchem_data(df: pd.DataFrame) -> pd.DataFrame:
    """Augment ChEMBL records with PubChem information.

    For each canonical SMILES string in ``df``, the function looks up the
    corresponding PubChem CID and basic chemical properties. The PubChem
    fields are appended to the input frame. If a SMILES string cannot be
    resolved, empty values are inserted.

    Parameters
    ----------
    df:
        Data frame returned by :func:`library.chembl_library.get_testitem`.

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
        cid = pl.get_cid_from_smiles(smi) or ""
        first_cid = cid.split("|")[0] if cid else ""
        if first_cid:
            props = pl.get_properties(first_cid)
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


def run_chembl(args: argparse.Namespace) -> int:
    """Execute compound retrieval from the ChEMBL API and augment with PubChem data.

    Parameters
    ----------
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
            sep=args.sep,
            encoding=args.encoding,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    logger.info("Retrieved %d identifiers", len(ids))
    logger.info("Fetching ChEMBL data in chunks of %d", args.chunk_size)
    df = cl.get_testitem(ids, chunk_size=args.chunk_size)
    logger.info("Retrieved %d rows from ChEMBL", len(df))
    logger.info("Augmenting results with PubChem data")
    df = add_pubchem_data(df)
    logger.info("PubChem augmentation completed")
    output = args.output_csv or io.default_output_path(args.input_csv)
    try:
        io.write_csv(df, output, sep=args.sep, encoding=args.encoding)
        logger.info("Wrote %d rows to %s", len(df), output)
        return 0
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="ChEMBL and PubChem compound data utilities"
    )
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    parser.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=Path("input.csv"),
        help="CSV file containing molecule identifiers",
    )
    parser.add_argument(
        "--output",
        dest="output_csv",
        type=Path,
        default=None,
        help="Destination CSV file (default: auto-generate)",
    )
    parser.add_argument(
        "--column",
        default="molecule_chembl_id",
        help="Column name in the input CSV containing identifiers",
    )
    parser.add_argument("--sep", default=",", help="CSV delimiter")
    parser.add_argument("--encoding", default="utf8", help="File encoding")
    parser.add_argument(
        "--chunk-size", type=int, default=5, help="Maximum number of IDs per request"
    )
    parser.set_defaults(func=run_chembl)
    return parser


def configure_logging(level: str) -> None:
    """Configure basic logging."""
    numeric_level = getattr(logging, level.upper(), logging.INFO)
    logging.basicConfig(level=numeric_level)


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
