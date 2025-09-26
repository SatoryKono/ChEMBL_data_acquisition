"""Command line interface for retrieving ChEMBL and PubChem compound data."""

from __future__ import annotations

import sys

# ruff: noqa: E402
from pathlib import Path

if __package__ is None:  # running as a script
    sys.path.append(str(Path(__file__).resolve().parents[1]))

import argparse
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from itertools import islice

import pandas as pd
import requests
from pandera.errors import SchemaErrors

from library import chembl_library as cl
from library import io
from library import molecule_catalog
from library import pubchem_library as pl
from library.molecule_catalog import load_parent_catalog
from library.chembl_client import ChemblClient
from library.cli import (
    LoggerConfig,
    apply_config_overrides,
    configure_logger,
)
from library.cli import (
    build_parser as base_parser,
)
from library.config import (
    ApiCfg,
    Config,
    MoleculeCatalogCfg,
    PubChemCfg,
    _serialize_paths,
    ensure_dirs,
    print_config,
)
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.pipeline_metadata import add_pipeline_metadata
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from library.validation import validate_testitems
from schemas import TestitemsSchema, normalize_testitems

_TYPO_PARENT_COLUMN = "parant_molecule_id"


def ensure_no_parant_column(df: pd.DataFrame) -> None:
    """Raise a :class:`ValueError` if the legacy typo column is present."""

    if _TYPO_PARENT_COLUMN in df.columns:
        raise ValueError(
            "unexpected column 'parant_molecule_id'; use 'parent_molecule_id' instead"
        )


PARENT_LOOKUP_SOURCE_CACHE = "cache"
PARENT_LOOKUP_SOURCE_REMOTE = "remote"
PARENT_LOOKUP_SOURCE_SKIPPED = "skipped"


@dataclass(frozen=True)
class ParentLookupStats:
    """Summary information about parent molecule enrichment."""

    source: str
    missing: int
    unique: int
    attached: int


def _cache_state(path: Path) -> tuple[bool, float | None]:
    """Return cache file presence and modification time."""

    try:
        stat = path.stat()
    except FileNotFoundError:
        return False, None
    return True, stat.st_mtime


def _resolve_parent_source(
    before: tuple[bool, float | None], after: tuple[bool, float | None]
) -> str:
    """Determine parent lookup source from cache state transitions."""

    before_exists, before_mtime = before
    after_exists, after_mtime = after

    if after_exists and before_exists and after_mtime == before_mtime:
        return PARENT_LOOKUP_SOURCE_CACHE
    if after_exists or before_exists != after_exists:
        return PARENT_LOOKUP_SOURCE_REMOTE
    return PARENT_LOOKUP_SOURCE_REMOTE


def _normalise_chembl_ids(series: pd.Series) -> pd.Series:
    """Return ``series`` normalised to upper-case ChEMBL identifiers."""

    normalised = (
        series.fillna("")
        .astype("string")
        .str.strip()
        .str.upper()
    )
    return normalised


def attach_parent_molecule_ids(
    df: pd.DataFrame,
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    catalog_cfg: MoleculeCatalogCfg,
    timeout: float | None,
    catalog: Mapping[str, str] | None = None,
    source: str | None = None,
) -> tuple[pd.DataFrame, ParentLookupStats]:
    """Attach parent molecule identifiers using the ChEMBL catalogue."""

    if "parant_molecule_id" in df.columns:
        raise ValueError(
            "Input frame contains unexpected column 'parant_molecule_id'."
        )

    result = df.copy()

    if result.empty:
        if catalog_cfg.parent_field not in result.columns:
            result[catalog_cfg.parent_field] = pd.Series(
                pd.NA, index=result.index, dtype="string"
            )
        stats = ParentLookupStats(
            source=PARENT_LOOKUP_SOURCE_SKIPPED,
            missing=0,
            unique=0,
            attached=0,
        )
        return result, stats

    child_column = catalog_cfg.child_field
    parent_column = catalog_cfg.parent_field

    if child_column not in result.columns:
        logger.warning("parent_lookup_missing_child_column", column=child_column)
        result[parent_column] = pd.Series(pd.NA, index=result.index, dtype="string")
        stats = ParentLookupStats(
            source=PARENT_LOOKUP_SOURCE_SKIPPED,
            missing=len(result),
            unique=0,
            attached=0,
        )
        return result, stats

    source_resolved = source
    catalog_data: dict[str, str]

    if catalog is None:
        cache_before = _cache_state(catalog_cfg.cache_path)
        catalog_data = load_parent_catalog(
            client=client,
            api_cfg=api_cfg,
            catalog_cfg=catalog_cfg,
            timeout=timeout,
        )
        cache_after = _cache_state(catalog_cfg.cache_path)
        source_resolved = _resolve_parent_source(cache_before, cache_after)
    else:
        catalog_data = dict(catalog)
        if source_resolved is None:
            cache_exists = catalog_cfg.cache_path.is_file()
            source_resolved = (
                PARENT_LOOKUP_SOURCE_CACHE
                if cache_exists
                else PARENT_LOOKUP_SOURCE_REMOTE
            )

    normalised_child = _normalise_chembl_ids(result[child_column])
    unique_children = normalised_child[normalised_child != ""].unique()

    parent_map = {
        key: catalog_data[key]
        for key in unique_children
        if key in catalog_data
    }
    missing_ids = [key for key in unique_children if key not in parent_map]
    fetched_remote = False

    if missing_ids and catalog is None:
        try:
            fetched = molecule_catalog.fetch_parent_catalog_for(
                missing_ids,
                client=client,
                api_cfg=api_cfg,
                timeout=timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.warning("parent_lookup_partial_fetch_failed", error=str(exc))
            fetched = {}
        else:
            if fetched:
                fetched_remote = True
        if fetched:
            catalog_data.update(fetched)
            parent_map.update(fetched)

    parent_series = normalised_child.map(parent_map).astype("string")

    if parent_column in result.columns:
        existing_parent = result[parent_column].astype("string")
    else:
        existing_parent = pd.Series(pd.NA, index=result.index, dtype="string")

    needs_fill = existing_parent.fillna("").str.strip() == ""
    combined_parent = existing_parent.copy()
    combined_parent.loc[needs_fill] = parent_series.loc[needs_fill]
    combined_parent = combined_parent.fillna(parent_series)
    result[parent_column] = combined_parent

    missing = int(combined_parent.isna().sum())
    attached = len(result) - missing

    stats = ParentLookupStats(
        source=(
            PARENT_LOOKUP_SOURCE_REMOTE
            if fetched_remote
            else source_resolved
        ),
        missing=missing,
        unique=int(len(unique_children)),
        attached=int(attached),
    )

    logger.info(
        "parent_lookup_progress",
        source=stats.source,
        unique=stats.unique,
        attached=stats.attached,
        missing=stats.missing,
    )

    return result, stats


def add_pubchem_data(df: pd.DataFrame, cfg: PubChemCfg) -> pd.DataFrame:
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
    if df.empty or "canonical_smiles" not in df.columns:
        return df

    smiles_list = df["canonical_smiles"].fillna("").tolist()
    # ``dict.fromkeys`` preserves the order of first occurrence while
    # removing duplicates. This allows progress output to reflect the
    # deterministic iteration order of SMILES strings.
    unique_smiles = [s for s in dict.fromkeys(smiles_list) if s]

    total = len(unique_smiles)
    if total:
        logger.info("pubchem_start", total=total)
    else:
        logger.info("pubchem_no_smiles")

    records: dict[str, dict[str, str]] = {}
    for idx, smi in enumerate(unique_smiles, start=1):
        logger.info("pubchem_progress", current=idx, total=total)
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
    limit = cfg.testitem.limit
    if limit is not None and limit < 0:
        logger.error(
            "invalid_limit",
            section="testitem.limit",
            limit=limit,
        )
        return 1

    # Initialise HTTP sessions for downstream HTTP calls
    pl.init_session(cfg.api, cfg.retry)
    # Initialise HTTP session for subsequent ChEMBL requests
    parent_stats = ParentLookupStats(
        source=PARENT_LOOKUP_SOURCE_SKIPPED,
        missing=0,
        unique=0,
        attached=0,
    )
    with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
        try:
            ids_iter = io.read_ids(
                args.input_csv, column=cfg.testitem.column, cfg=cfg.io
            )
            if limit is not None:
                ids = list(islice(ids_iter, limit))
                logger.info("process_limit", limit=len(ids))
            else:
                ids = list(ids_iter)
        except (FileNotFoundError, ValueError) as exc:
            logger.error(
                "read_fail",
                error=str(exc),
                path=str(args.input_csv),
            )
            return 1

        logger.info("identifiers_retrieved", count=len(ids))
        logger.info("chembl_fetch_start", batch_size=cfg.testitem.batch_size)

        cache_before = _cache_state(cfg.molecule_catalog.cache_path)

        try:
            parent_catalog = load_parent_catalog(
                client=client,
                api_cfg=cfg.api,
                catalog_cfg=cfg.molecule_catalog,
                timeout=cfg.testitem.timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error(
                "parent_catalog_invalid",
                error=str(exc),
                path=str(cfg.molecule_catalog.cache_path),
            )
            return 1

        cache_after = _cache_state(cfg.molecule_catalog.cache_path)
        parent_catalog_source = _resolve_parent_source(cache_before, cache_after)

        try:
            df = cl.get_testitem(
                ids,
                cfg=cfg.api,
                client=client,
                chunk_size=cfg.testitem.batch_size,
                timeout=cfg.testitem.timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error(
                "testitem_fetch_failed",
                error=str(exc),
                batch_size=cfg.testitem.batch_size,
                timeout=cfg.testitem.timeout,
            )
            return 1
        logger.info("chembl_fetch_done", rows=len(df))
        parent_column = cfg.molecule_catalog.parent_field
        if parent_catalog and "molecule_chembl_id" in df.columns:
            normalised_ids = _normalise_chembl_ids(df["molecule_chembl_id"])
            mapped = normalised_ids.map(parent_catalog)
            if parent_column in df.columns:
                df[parent_column] = df[parent_column].fillna(mapped)
            else:
                df[parent_column] = mapped
        logger.info("pubchem_augment_start")
        df = add_pubchem_data(df, cfg.pubchem)
        logger.info("pubchem_augment_done")
 
        try:
            ensure_no_parant_column(df)
        except ValueError as exc:
            logger.error(
                "invalid_column",
                column=_TYPO_PARENT_COLUMN,
                error=str(exc),
            )
            return 1
 
        logger.info("parent_lookup_start")
        try:
            df, parent_stats = attach_parent_molecule_ids(
                df,
                client=client,
                api_cfg=cfg.api,
                catalog_cfg=cfg.molecule_catalog,
                timeout=cfg.testitem.timeout,
                catalog=parent_catalog,
                source=parent_catalog_source,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error("parent_lookup_failed", error=str(exc))
            return 1
        logger.info(
            "parent_lookup_done",
            source=parent_stats.source,
            unique=parent_stats.unique,
            attached=parent_stats.attached,
            missing=parent_stats.missing,
        )
 
        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        df = normalize_testitems(df)
        df = add_pipeline_metadata(df)
        # Determine column order: schema columns first, followed by
        # additional fields sorted alphabetically.
        schema_cols = list(TestitemsSchema.columns)
        head = [c for c in schema_cols if c in df.columns]
        tail = sorted(c for c in df.columns if c not in schema_cols)
        col_order = head + tail
        rows_total = len(df)
        exit_code = 0
    required_cols = {
        name for name, col in TestitemsSchema.columns.items() if col.required
    }
    optional_cols = set(TestitemsSchema.columns) - required_cols
    missing_required = required_cols - set(df.columns)
    missing_optional = optional_cols - set(df.columns)
    if not missing_required:
        if missing_optional:
            logger.warning(
                "optional_columns_missing",
                columns=sorted(missing_optional),
            )
        try:
            validation_result = validate_testitems(df, return_result=True)
        except SchemaErrors as exc:
            failure_path = Path(output).with_name(
                f"{Path(output).stem}_failure_cases.csv"
            )
            errors = SidecarErrors()
            for row in exc.failure_cases.to_dict("records"):
                errors.add_error(row)
            errors.save(failure_path)
            logger.error(
                "validation_failed",
                failures=len(exc.failure_cases),
                path=str(failure_path),
            )
            df = getattr(exc, "validated_data", df)
            exit_code = 1
        else:
            df = validation_result.data
            if not validation_result.failure_cases.empty:
                failure_path = Path(output).with_name(
                    f"{Path(output).stem}_failure_cases.csv"
                )
                errors = SidecarErrors()
                for row in validation_result.failure_cases.to_dict("records"):
                    errors.add_error(row)
                errors.save(failure_path)
                logger.error(
                    "validation_failed",
                    failures=len(validation_result.failure_cases),
                    path=str(failure_path),
                )
                exit_code = 1
    else:
        logger.warning(
            "validation_skipped",
            missing_columns=sorted(missing_required),
        )
    rows_kept = len(df)
    rows_dropped = rows_total - rows_kept
    try:
        key_cols = ["molecule_chembl_id"]
        csv_path = io.write_csv(
            df,
            output,
            cfg=cfg,
            key_cols=key_cols or None,
            col_order=col_order,
        )
        logger.info("write_done", rows=rows_kept, path=str(csv_path))
    except OSError as exc:
        logger.error(
            "write_fail",
            error=str(exc),
            path=str(output),
        )
        return 1

    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_kept,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(csv_path),
        "parent_lookup_source": parent_stats.source,
        "parent_lookup_missing": parent_stats.missing,
    }
    write_meta_yaml(
        csv_path=csv_path,
        command=" ".join(sys.argv),
        config_subset=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(args.input_csv)},
        stats=stats,
        schema="TestitemsSchema",
    )
    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error(
            "quality_report_failed",
            error=str(exc),
            path=str(output),
        )
        return 1
    return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser, log_cfg = base_parser(
        "ChEMBL and PubChem compound data utilities",
        column="molecule_chembl_id",
        chunk_size=5,
        size_option="--batch-size",
        size_dest="batch_size",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be a positive integer")
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "timeout": "testitem.timeout",
                "column": "testitem.column",
                "batch_size": "testitem.batch_size",
                "limit": "testitem.limit",
            },
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error(
            "directory_setup_failed",
            error=str(exc),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    exit_code: int = args.func(cfg, args)
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
