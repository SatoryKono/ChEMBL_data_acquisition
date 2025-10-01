"""Command line interface for retrieving document metadata from external sources.

The tool integrates :mod:`library.pubmed_library` and
:mod:`library.chembl_library` to collect information about publications from
several public APIs.  The interface mirrors :mod:`scripts.get_target_data` and
provides three sub-commands:

``pubmed``
    Query PubMed, Semantic Scholar, OpenAlex and CrossRef for a list of PMIDs.
``chembl``
    Retrieve document information from the ChEMBL API.
``all``
    Run the ChEMBL and PubMed pipelines and merge the results.

Example
-------
Fetch PubMed metadata for identifiers listed in ``pmids.csv``::

    python scripts/get_document_data.py pubmed --config config/config.yaml --input pmids.csv --output output.csv

The input file must contain a ``PMID`` column.

"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

import requests


def _ensure_project_root() -> None:
    """Ensure the repository root is discoverable when executed as a script."""

    script_path = Path(__file__).resolve()
    project_root = script_path.parents[1]
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)


if __package__ in {None, ""}:
    _ensure_project_root()

from library import cli
from library import document_postprocessing as dp
from library import io
from library import openalex_crossref_library as ocl
from library import pubmed_library as pl
from library import semantic_scholar_library as ssl
import library.document_pipeline as document_pipeline
from library.clients import ChemblClient, _chunked
from library.cli import (
    LoggerConfig,
    build_root_parser,
    configure_logger,
    path_argument,
    positive_int,
)
from library.config import Config, ensure_dirs, print_config, session_with_retry
from library.document_pipeline import (
    DOCUMENT_EXPORT_COLUMNS,
    EXPORT_COALESCE_SOURCES,
    EXPORT_COLUMN_RENAMES,
    EXPORT_SORT_FALLBACK,
    EXPORT_STREAM_CHUNK_SIZE,
    NUMERIC_EXPORT_COLUMNS,
    coalesce_columns,
    build_fallback_doi_map,
    fetch_pubmed_records,
    finalise_document_export,
    iter_document_export_chunks,
    load_export_ready_frame,
    prepare_document_export_frame,
    resolve_export_chunk_size,
    run_all_pipeline,
    run_chembl_pipeline,
    run_pubmed_pipeline,
    write_document_export_chunks,
)
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.pipeline_metadata import add_pipeline_metadata
from library.rate_limiter import get_limiter
from library.table_quality import analyze_table_quality
from library.log import logger
from schemas import DocumentsSchema, normalize_documents


# Backwards-compatible re-exports for unit tests targeting script internals.
_build_fallback_doi_map = build_fallback_doi_map
_fetch_pubmed_records = fetch_pubmed_records
_finalise_export = finalise_document_export
_iter_export_chunks = iter_document_export_chunks
_coalesce_columns = coalesce_columns
_load_export_ready_frame = load_export_ready_frame
_prepare_export_frame = prepare_document_export_frame
_resolve_chunk_size = resolve_export_chunk_size
_write_export_chunks = write_document_export_chunks
_NUMERIC_EXPORT_COLUMNS = NUMERIC_EXPORT_COLUMNS
_EXPORT_COLUMNS = DOCUMENT_EXPORT_COLUMNS
_EXPORT_COLUMN_RENAMES = EXPORT_COLUMN_RENAMES
_EXPORT_COALESCE_SOURCES = EXPORT_COALESCE_SOURCES
_EXPORT_SORT_FALLBACK = EXPORT_SORT_FALLBACK
_EXPORT_STREAM_CHUNK_SIZE = EXPORT_STREAM_CHUNK_SIZE
ThreadPoolExecutor = document_pipeline.ThreadPoolExecutor
as_completed = document_pipeline.as_completed


def run_pubmed(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the PubMed pipeline using library helpers."""

    return run_pubmed_pipeline(cfg, args)


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ChEMBL pipeline using library helpers."""

    return run_chembl_pipeline(cfg, args)


def run_all(cfg: Config, args: argparse.Namespace) -> int:
    """Run both pipelines and merge their outputs."""

    return run_all_pipeline(cfg, args)


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the argument parser for document utilities."""

    root, shared, log_cfg = build_root_parser()
    parser = argparse.ArgumentParser(
        description="Document data utilities", parents=[root]
    )
    sub = parser.add_subparsers(dest="command", required=True)

    pubmed = sub.add_parser(
        "pubmed", parents=[shared], help="Fetch data from PubMed and related APIs"
    )
    pubmed.add_argument(
        "--column", default="PMID", help="Column name containing identifiers"
    )
    pubmed.add_argument(
        "--sleep", type=float, default=5.0, help="Seconds to sleep between requests"
    )
    pubmed.add_argument(
        "--workers", type=int, default=1, help="Number of concurrent requests"
    )
    pubmed.add_argument(
        "--batch-size",
        type=positive_int,
        default=100,
        help="Maximum PMIDs per PubMed request",
    )
    pubmed.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    pubmed.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    pubmed.add_argument(
        "--openalex-rps",
        type=float,
        default=None,
        help="Requests per second limit for OpenAlex",
    )
    pubmed.add_argument(
        "--crossref-rps",
        type=float,
        default=None,
        help="Requests per second limit for CrossRef",
    )
    pubmed.add_argument(
        "--fallback-doi-csv",
        type=path_argument,
        default=None,
        help="Optional CSV file providing PMID to DOI overrides",
    )
    pubmed.add_argument(
        "--fallback-doi-pmid-column",
        default="PMID",
        help="Column containing PubMed identifiers in fallback CSV",
    )
    pubmed.add_argument(
        "--fallback-doi-value-column",
        default="DOI",
        help="Column containing DOI values in fallback CSV",
    )
    pubmed.set_defaults(func=run_pubmed)

    chembl = sub.add_parser(
        "chembl", parents=[shared], help="Fetch document information from ChEMBL"
    )
    chembl.add_argument(
        "--column",
        default="document_chembl_id",
        help="Column name containing identifiers",
    )
    chembl.add_argument(
        "--chunk-size",
        type=positive_int,
        default=5,
        help="Maximum number of IDs per request",
    )
    chembl.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    chembl.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    chembl.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    chembl.set_defaults(func=run_chembl)

    all_cmd = sub.add_parser(
        "all", parents=[shared], help="Run both ChEMBL and PubMed pipelines"
    )
    all_cmd.add_argument(
        "--column",
        default="document_chembl_id",
        help="Column in the input CSV",
    )
    all_cmd.add_argument(
        "--chunk-size",
        type=positive_int,
        default=5,
        help="Maximum IDs per request",
    )
    all_cmd.add_argument(
        "--sleep",
        type=float,
        default=5.0,
        help="Seconds to sleep between PubMed requests",
    )
    all_cmd.add_argument(
        "--workers", type=int, default=1, help="Number of concurrent PubMed requests"
    )
    all_cmd.add_argument(
        "--batch-size",
        type=positive_int,
        default=50,
        help="Maximum PMIDs per PubMed request",
    )
    all_cmd.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    all_cmd.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    all_cmd.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    all_cmd.add_argument(
        "--openalex-rps",
        type=float,
        default=None,
        help="Requests per second limit for OpenAlex",
    )
    all_cmd.add_argument(
        "--crossref-rps",
        type=float,
        default=None,
        help="Requests per second limit for CrossRef",
    )
    all_cmd.set_defaults(func=run_all)

    parser.subparsers_map = {  # type: ignore[attr-defined]
        "pubmed": pubmed,
        "chembl": chembl,
        "all": all_cmd,
    }

    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""

    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)

    subparser_map = getattr(parser, "subparsers_map", {})
    subparser = subparser_map.get(args.command, parser)
    limit_value = getattr(args, "limit", None)
    if limit_value is not None and limit_value <= 0:
        subparser.error("--limit must be a positive integer")
    offset_value = getattr(args, "offset", 0)
    if offset_value < 0:
        subparser.error("--offset must be zero or a positive integer")

    log_cfg.level = args.log_level
    logger_instance = configure_logger(log_cfg)
    logger_instance.info("pipeline_start", run_id=log_cfg.run_id)

    mapping = {
        "column": f"document.{args.command}.column",
        "limit": f"document.{args.command}.limit",
    }
    if args.command == "pubmed":
        mapping.update(
            {
                "sleep": "document.pubmed.sleep",
                "workers": "document.pubmed.workers",
                "batch_size": "document.pubmed.batch_size",
            }
        )
    elif args.command == "chembl":
        mapping.update(
            {
                "chunk_size": "document.chembl.chunk_size",
                "timeout": "document.chembl.timeout",
            }
        )
    elif args.command == "all":
        mapping.update(
            {
                "chunk_size": "document.all.chunk_size",
                "sleep": "document.all.sleep",
                "workers": "document.all.workers",
                "batch_size": "document.all.batch_size",
                "timeout": "document.all.timeout",
            }
        )
    try:
        cfg: Config = cli.apply_config_overrides(
            args,
            subparser,
            args.config,
            mapping=mapping
            | {
                "openalex_rps": "openalex.rps",
                "crossref_rps": "crossref.rps",
            },
            base_parser=parser,
        )
        cfg.api.timeout_read = getattr(args, "timeout", cfg.api.timeout_read)
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger_instance = configure_logger(log_cfg)
    except (ValueError, TypeError) as exc:
        logger.error("config_override_error", error=str(exc))
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("directory_setup_failed", error=str(exc))
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    exit_code = getattr(args, "func", run_pubmed)(cfg, args)
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
