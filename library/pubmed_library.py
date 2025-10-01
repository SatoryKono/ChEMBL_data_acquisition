"""Utilities for retrieving and merging publication metadata.

This module exposes a command line interface and re-exports the
implementation located in :mod:`library.pubmed`.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from datetime import date
import sys
from pathlib import Path

import pandas as pd

from . import cli
from .clients.semantic_scholar import (
    fetch_semantic_scholar,
    fetch_semantic_scholar_batch,
)
from .cli import LoggerConfig, configure_logger, path_argument
from .cli import build_parser as base_parser
from .config import (
    Config,
    _serialize_paths,
    ensure_dirs,
    print_config,
    session_with_retry,
)
from .csv_utils import write_csv_deterministic
from .clients.pubmed import PubMedClient
from .pubmed import (
    EMPTY_PUBMED,
    combine,
    fetch_crossref,
    fetch_openalex,
    fetch_pubmed,
    fetch_pubmed_batch,
    find_all,
    find_one,
    merge_records,
    parse_pubmed_article,
    print_results,
    read_pmids,
    text_or_none,
)
from .rate_limiter import get_limiter
from .metadata import Stats, file_sha256, write_meta_yaml as write_pipeline_meta_yaml
from .sidecar import SidecarErrors
from .table_quality import analyze_table_quality

__all__ = [
    "Config",
    "read_pmids",
    "fetch_pubmed_batch",
    "fetch_pubmed",
    "fetch_semantic_scholar",
    "fetch_semantic_scholar_batch",
    "fetch_openalex",
    "fetch_crossref",
    "text_or_none",
    "combine",
    "find_one",
    "find_all",
    "parse_pubmed_article",
    "EMPTY_PUBMED",
    "merge_records",
    "print_results",
    "parse_args",
    "main",
]


REQUIRED_FIELDS = ("PubMed.PMID",)


def parse_args(
    argv: Sequence[str] | None = None,
) -> tuple[argparse.Namespace, argparse.ArgumentParser, LoggerConfig]:
    """Parse command-line arguments."""
    parser, log_cfg = base_parser("Fetch publication metadata by PMID", column="PMID")
    parser.add_argument(
        "--input-csv",
        dest="input_csv",
        type=path_argument,
        default=argparse.SUPPRESS,
        help="Input CSV path with PMID column",
    )
    parser.add_argument(
        "--keep-verbose-dumps",
        action="store_true",
        help="Log combined metadata dumps at INFO level for troubleshooting",
    )
    args = parser.parse_args(argv)
    return args, parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line interface entry point."""
    args, parser, log_cfg = parse_args(argv)
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg = cli.apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "column": "document.pubmed.column",
                "chunk_size": "document.pubmed.batch_size",
            },
        )
        cfg = Config.model_validate(cfg.model_dump())
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg)
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

    pubmed_rps = cfg.pubmed.rps or cfg.rate.global_rps
    pubmed_burst = cfg.pubmed.burst or cfg.rate.global_burst
    limiter = get_limiter("pubmed", pubmed_rps, pubmed_burst)
    delay = 1.0 / pubmed_rps if pubmed_rps > 0 else 0.0

    raw_argv = list(argv) if argv is not None else sys.argv[1:]
    command = " ".join(["library.pubmed_library", *map(str, raw_argv)])
    pmid_df = read_pmids(args.input_csv, cfg=cfg.pubmed)
    pmids = pmid_df["PMID"].tolist()
    openalex_limiter = get_limiter("openalex", cfg.openalex.rps, cfg.openalex.burst)
    crossref_limiter = get_limiter("crossref", cfg.crossref.rps, cfg.crossref.burst)
    pubmed_client = PubMedClient(cfg.pubmed)

    valid_records: list[dict[str, str]] = []
    errors = SidecarErrors()
    rows_total = 0
    rows_kept = 0
    total_failures = 0
    batch_size = cfg.document.pubmed.batch_size
    dump_level = "INFO" if args.keep_verbose_dumps else "DEBUG"
    with session_with_retry(cfg.api, cfg.retry) as session:
        for i in range(0, len(pmids), batch_size):
            batch_pmids = pmids[i : i + batch_size]
            limiter.acquire()
            pubmed_list = fetch_pubmed_batch(
                session, batch_pmids, delay, client=pubmed_client
            )
            limiter.acquire()
            semsch_list = fetch_semantic_scholar_batch(
                session, batch_pmids, delay, cfg=cfg.semantic_scholar
            )
            semsch_map = {s.get("scholar.PMID"): s for s in semsch_list}
            for index, pubmed in enumerate(pubmed_list):
                source_pmid = batch_pmids[index]
                pmid = pubmed.get("PubMed.PMID", "")
                resolved_pmid = pmid or source_pmid
                semsch = semsch_map.get(pmid) or semsch_map.get(source_pmid, {})

                openalex = fetch_openalex(
                    session,
                    resolved_pmid,
                    cfg=cfg.openalex,
                    limiter=openalex_limiter,
                )
                doi = pubmed.get("PubMed.DOI") or semsch.get("scholar.DOI") or ""
                crossref = fetch_crossref(
                    session, doi, cfg=cfg.crossref, limiter=crossref_limiter
                )

                combined = merge_records(pubmed, semsch, openalex, crossref)
                print_results([combined], level=dump_level)
                rows_total += 1
                missing_fields = [
                    field
                    for field in REQUIRED_FIELDS
                    if not str(combined.get(field, "")).strip()
                ]
                if missing_fields:
                    total_failures += 1
                    errors.add_error(
                        {
                            "row_index": rows_total - 1,
                            "missing_fields": "|".join(missing_fields),
                            "source_pmid": source_pmid,
                            "PubMed.PMID": combined.get("PubMed.PMID", ""),
                        }
                    )
                    continue
                valid_records.append(combined)
                rows_kept += 1

    df = pd.DataFrame.from_records(valid_records)
    output_path = (
        Path(args.output_csv)
        if args.output_csv
        else Path(f"output.{Path(args.input_csv).stem}_{date.today():%Y%m%d}.csv")
    )
    failure_path = output_path.with_name(f"{output_path.stem}_failure_cases.csv")
    csv_path = write_csv_deterministic(df, output_path, key_cols=sorted(df.columns))
    rows_dropped = rows_total - rows_kept
    if total_failures:
        logger.error(
            "validation_failed",
            failures=total_failures,
            path=str(failure_path),
        )
        errors.save(failure_path)
    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_kept,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(csv_path),
    }
    write_pipeline_meta_yaml(
        csv_path=csv_path,
        command=command,
        config_subset=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(args.input_csv)},
        stats=stats,
        schema="PubMedMetadata",
    )
    try:
        analyze_table_quality(df, table_name=str(csv_path.with_suffix("")))
    except ValueError as exc:
        logger.error(
            "quality_report_failed",
            error=str(exc),
            path=str(csv_path),
        )
        return 1
    logger.info("file_written", path=str(csv_path))
    logger.info("pipeline_done", run_id=log_cfg.run_id)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
