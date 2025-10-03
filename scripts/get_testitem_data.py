"""Command line interface for retrieving ChEMBL test item data.

The module wraps :func:`library.testitem_pipeline.run_testitem_pipeline` while
exposing helpers that tests can import directly. Entry points return numeric
exit codes rather than terminating the interpreter to simplify orchestration.
The :func:`ensure_no_parant_column` helper guards against legacy CSV exports
that still include the misspelled parent identifier column.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Hashable, MutableMapping, Sequence, cast

import pandas as pd


# ===== Parameters =====

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_NAME = "testitem.csv"
DEFAULT_OUTPUT_STEM = "testitems"


if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from library.utils.bootstrap import ensure_project_root


if __package__ in {None, ""}:
    ensure_project_root()

from library import cli  # noqa: F401 - re-exported for monkeypatching in tests
from library import io
from library.integration import molecule_catalog
from library.integration import pubchem_library as pl
from library.cli import LoggerConfig
from library.cli import build_parser as base_parser
from library.cli_utils import run_cli_command
from library.config import (
    ApiCfg,
    Config,
    IoCfg,
    MoleculeCatalogCfg,
    PubChemCfg,
)
from library.common.log import logger
from library.clients import pubchem as pc  # noqa: F401 - patched in tests
import library.testitem_pipeline as pipeline
from library.integration.chembl_client import ChemblClient
from library.testitem_pipeline import (
    PUBCHEM_CID_CACHE_ENCODING,
    PUBCHEM_COLUMNS,
    ReadInputIdsResult,
    TestitemPipelineOptions,
    _DEFAULT_CATALOG_CFG,
    _FETCH_ERROR_SAMPLE_SIZE,
    _MOLECULE_HIERARCHY_COLUMNS,
    _PUBCHEM_CACHE_SCHEMA_VERSION,
    _TYPO_PARENT_COLUMN,
    analyze_table_quality,
    ensure_no_parant_column,
    file_sha256,
    fetch_testitems,
    integrate_missing_identifiers,
    load_parent_catalog,
    query_parent_catalog,
    read_input_ids,
    run_testitem_pipeline,
    update_parent_catalog_cache,
    write_meta_yaml,
    write_parent_catalog_cache,
    _prepare_pubchem_api_cfg,
    _write_pubchem_cid_cache,
    PARENT_LOOKUP_SOURCE_CACHE,
    PARENT_LOOKUP_SOURCE_LOOKUP,
    PARENT_LOOKUP_SOURCE_PARTIAL,
    PARENT_LOOKUP_SOURCE_SKIPPED,
    PARENT_LOOKUP_SOURCE_SYNC,
)
from library.testitem_pipeline import catalog as pipeline_catalog
from library.testitem_pipeline import pubchem as pipeline_pubchem

LoadMoleculeHierarchyLookup = pipeline.LoadMoleculeHierarchyLookup
load_molecule_hierarchy_lookup = pipeline.load_molecule_hierarchy_lookup
attach_parent_molecule_ids = pipeline.attach_parent_molecule_ids

_normalise_identifier = pipeline_pubchem._normalise_identifier
_pubchem_identifiers = pipeline_pubchem._pubchem_identifiers
_pubchem_resolution_key = pipeline_pubchem._pubchem_resolution_key
_load_pubchem_cid_cache = pipeline_pubchem._load_pubchem_cid_cache
resolve_pubchem_cid = pipeline_pubchem.resolve_pubchem_cid


_normalise_parent_identifier = pipeline_catalog._normalise_parent_identifier
_load_molecule_hierarchy_mapping = pipeline_catalog._load_molecule_hierarchy_mapping

def add_pubchem_data(
    df: pd.DataFrame,
    cfg: PubChemCfg,
    *,
    client: ChemblClient | None = None,
    api_cfg: ApiCfg | None = None,
    timeout: float | None = None,
    cid_cache: MutableMapping[str, str | None] | None = None,
    resolution_cache: MutableMapping[Hashable, pl.PubChemResolution] | None = None,
    parent_record_cache: MutableMapping[str, pd.Series | None] | None = None,
    testitem_fields: Sequence[str] | None = None,
    request_limit: int = 1000,
) -> pd.DataFrame:
    """Augment ChEMBL records with PubChem information.

    Delegates to :func:`library.testitem_pipeline.add_pubchem_data` while
    relaxing the ``resolution_cache`` type to align with
    :func:`library.integration.pubchem_library.resolve_pubchem_record`.
    """

    return pipeline.add_pubchem_data(
        df,
        cfg,
        client=client,
        api_cfg=api_cfg,
        timeout=timeout,
        cid_cache=cid_cache,
        resolution_cache=cast(
            MutableMapping[tuple[str | None, ...], pl.PubChemResolution] | None,
            resolution_cache,
        ),
        parent_record_cache=parent_record_cache,
        testitem_fields=testitem_fields,
        request_limit=request_limit,
    )


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Invoke the reusable test item pipeline with CLI parameters.

    Parameters
    ----------
    cfg : Config
        Application configuration containing ChEMBL, PubChem and IO settings.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate that identifier loading,
        network requests or CSV export failed inside the test item pipeline.
    """

    output_csv = getattr(args, "output_csv", None)
    options = TestitemPipelineOptions(
        input_csv=Path(args.input_csv),
        output_csv=Path(output_csv) if output_csv else None,
    )
    return run_testitem_pipeline(cfg, options)


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the test item pipeline handling ``--skip-existing`` semantics."""

    output_path = Path(
        args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    )
    args.output_csv = output_path
    if args.skip_existing and output_path.exists() and not args.force:
        logger.info("pipeline_skip_existing", output=str(output_path))
        return 0
    logger.info(
        "testitem_pipeline_start",
        input=str(args.input_csv),
        output=str(output_path),
        limit=getattr(cfg.testitem, "limit", None),
        offset=getattr(args, "offset", getattr(cfg.testitem, "offset", None)),
        batch_size=getattr(cfg.testitem, "batch_size", None),
        timeout=getattr(cfg.testitem, "timeout", None),
    )
    exit_code = run_chembl(cfg, args)
    if exit_code == 0:
        logger.info("testitem_pipeline_done", output=str(output_path))
    else:
        logger.error(
            "testitem_pipeline_failed",
            output=str(output_path),
            exit_code=exit_code,
        )
    return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        Parser configured with the common CLI options and the associated logging
        configuration used by :func:`main`.
    """

    parser, log_cfg = base_parser(
        "ChEMBL and PubChem compound data utilities",
        column="molecule_chembl_id",
        chunk_size=1000,
        size_option="--batch-size",
        size_dest="batch_size",
    )
    parser.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
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
        help=(
            "Maximum number of identifiers to process; use 0 to skip processing"
        ),
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults.

    Parameters
    ----------
    argv : Sequence[str] | None, optional
        Command-line arguments to parse. When ``None`` the values from
        :data:`sys.argv` are used.

    Returns
    -------
    int
        ``0`` when the pipeline finishes successfully, non-zero otherwise.

    Raises
    ------
    SystemExit
        Raised when invalid command-line options are provided.
    """

    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)

    cli.prepare_io_paths(
        args,
        input_default=DEFAULT_INPUT_NAME,
        output_stem=DEFAULT_OUTPUT_STEM,
    )

    if args.limit == 0:
        logger.info("pipeline_skip_limit", limit=args.limit)
        return 0

    if args.limit is not None and args.limit < 0:
        parser.error("--limit must be zero or a positive integer")
    if args.offset < 0:
        parser.error("--offset must be zero or a positive integer")

    mapping = {
        "timeout": "testitem.timeout",
        "column": "testitem.column",
        "batch_size": "testitem.batch_size",
        "limit": "testitem.limit",
        "offset": "testitem.offset",
    }
    return run_cli_command(
        args=args,
        parser=parser,
        log_cfg=log_cfg,
        mapping=mapping,
        run=run,
        logger=logger,
    )


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
