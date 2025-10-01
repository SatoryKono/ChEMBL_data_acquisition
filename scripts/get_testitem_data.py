"""Command line interface for retrieving ChEMBL and PubChem compound data."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence


_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from library import cli  # noqa: F401 - re-exported for monkeypatching in tests
from library.cli import LoggerConfig
from library.cli import build_parser as base_parser
from library.cli_utils import run_cli_command
from library.config import Config
from library.log import logger
from library.clients import pubchem as pc  # noqa: F401 - patched in tests
from library.testitem_pipeline import (
    ReadInputIdsResult as _ReadInputIdsResult,
    TestitemPipelineOptions,
    analyze_table_quality as _analyze_table_quality,
    fetch_testitems as _fetch_testitems,
    load_parent_catalog as _load_parent_catalog,
    query_parent_catalog as _query_parent_catalog,
    read_input_ids as _read_input_ids,
    run_testitem_pipeline,
)
from library.testitem_pipeline import (
    _prepare_pubchem_api_cfg as _pipeline_prepare_pubchem_api_cfg,
    _load_pubchem_cid_cache as _pipeline_load_pubchem_cid_cache,
    integrate_missing_identifiers as _integrate_missing_identifiers,
)

# Re-export helpers consumed directly by unit tests.
ReadInputIdsResult = _ReadInputIdsResult
read_input_ids = _read_input_ids
fetch_testitems = _fetch_testitems
load_parent_catalog = _load_parent_catalog
query_parent_catalog = _query_parent_catalog
_prepare_pubchem_api_cfg = _pipeline_prepare_pubchem_api_cfg
_load_pubchem_cid_cache = _pipeline_load_pubchem_cid_cache
analyze_table_quality = _analyze_table_quality
_integrate_missing_identifiers = _integrate_missing_identifiers


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Invoke the reusable test item pipeline with CLI parameters."""

    output_csv = getattr(args, "output_csv", None)
    options = TestitemPipelineOptions(
        input_csv=Path(args.input_csv),
        output_csv=Path(output_csv) if output_csv else None,
    )
    return run_testitem_pipeline(cfg, options)


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""

    parser, log_cfg = base_parser(
        "ChEMBL and PubChem compound data utilities",
        column="molecule_chembl_id",
        chunk_size=1000,
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
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""

    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)

    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be a positive integer")
    if args.offset < 0:
        parser.error("--offset must be zero or a positive integer")

    run_callable = getattr(args, "func", run_chembl)
    return run_cli_command(
        args=args,
        parser=parser,
        log_cfg=log_cfg,
        mapping={
            "timeout": "testitem.timeout",
            "column": "testitem.column",
            "batch_size": "testitem.batch_size",
            "limit": "testitem.limit",
            "offset": "testitem.offset",
        },
        run=run_callable,
        logger=logger,
    )


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
