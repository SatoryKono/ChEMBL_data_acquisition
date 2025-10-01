"""Command line interface for retrieving ChEMBL and PubChem compound data."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from library import cli
from library import molecule_catalog, pubchem_library as pl
from library import testitem_pipeline as pipeline
from library.cli import LoggerConfig
from library.cli import build_parser as base_parser
from library.cli.runner import run_cli_command
from library.config import Config
from library.log import logger

_CONFIG_OVERRIDES: dict[str, str] = {
    "timeout": "testitem.timeout",
    "column": "testitem.column",
    "batch_size": "testitem.batch_size",
    "limit": "testitem.limit",
    "offset": "testitem.offset",
}

load_parent_catalog = pipeline.load_parent_catalog
query_parent_catalog = pipeline.query_parent_catalog
update_parent_catalog_cache = pipeline.update_parent_catalog_cache
write_parent_catalog_cache = pipeline.write_parent_catalog_cache
analyze_table_quality = pipeline.analyze_table_quality
add_pubchem_data = pipeline.add_pubchem_data
load_molecule_hierarchy_lookup = pipeline.load_molecule_hierarchy_lookup
file_sha256 = pipeline.file_sha256
write_meta_yaml = pipeline.write_meta_yaml
read_input_ids = pipeline.read_input_ids
fetch_testitems = pipeline.fetch_testitems
prepare_parent_enrichment = pipeline.prepare_parent_enrichment
run_parent_enrichment = pipeline.run_parent_enrichment
finalize_output = pipeline.finalize_output


def _sync_pipeline_overrides() -> None:
    """Propagate script-level monkeypatches to the shared pipeline module."""

    pipeline.load_parent_catalog = load_parent_catalog
    pipeline.query_parent_catalog = query_parent_catalog
    pipeline.update_parent_catalog_cache = update_parent_catalog_cache
    pipeline.write_parent_catalog_cache = write_parent_catalog_cache
    pipeline.analyze_table_quality = analyze_table_quality
    pipeline.add_pubchem_data = add_pubchem_data
    pipeline.load_molecule_hierarchy_lookup = load_molecule_hierarchy_lookup
    pipeline.molecule_catalog.fetch_parent_catalog_for = (
        molecule_catalog.fetch_parent_catalog_for
    )
    pipeline.file_sha256 = file_sha256
    pipeline.write_meta_yaml = write_meta_yaml
    pipeline.read_input_ids = read_input_ids
    pipeline.fetch_testitems = fetch_testitems
    pipeline.prepare_parent_enrichment = prepare_parent_enrichment
    pipeline.run_parent_enrichment = run_parent_enrichment
    pipeline.finalize_output = finalize_output


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
    parser.set_defaults(func=_run_command)
    return parser, log_cfg


def _run_command(cfg: Config, args: argparse.Namespace) -> int:
    """Bridge parsed CLI arguments to the pipeline implementation."""

    _sync_pipeline_overrides()
    return pipeline.run_pipeline(
        cfg,
        input_csv=args.input_csv,
        output_csv=args.output_csv,
        offset=args.offset,
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""

    return run_cli_command(build_parser, argv, _CONFIG_OVERRIDES)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
