"""Utilities for retrieving and merging publication metadata.

This module exposes a command line interface and re-exports the
implementation located in :mod:`library.pubmed`.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Iterator, Sequence
from copy import deepcopy
from datetime import date
from pathlib import Path
from typing import cast

import pandas as pd

from library.schemas import DocumentsSchema

from .. import cli
from ..cli import LoggerConfig, configure_logger, path_argument
from ..cli import build_parser as base_parser
from ..cli.pipeline_definition import MetadataHook, PipelineDefinition, Validator
from ..cli_utils import run_pipeline
from ..clients.pubmed import PubMedClient
from ..clients.semantic_scholar import (
    fetch_semantic_scholar,
    fetch_semantic_scholar_batch,
)
from ..common.csv_utils import write_csv_chunks_deterministic
from ..common.rate_limiter import RateLimiter, get_limiter
from ..config.loader import _serialize_paths, ensure_dirs, print_config
from ..config.models import (
    ApiCfg,
    Config,
    ConfigError,
    CrossRefCfg,
    OpenAlexCfg,
    RetryCfg,
    SemanticScholarCfg,
)
from ..config.runtime import session_with_retry
from ..metadata import Stats, file_sha256, write_meta_yaml
from ..normalization import normalize_documents
from ..pipelines.common import add_pipeline_metadata
from ..pubmed import (
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
from ..qa.reporting import build_table_quality_hook
from ..validation import ValidationResult as SchemaValidationResult

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


_PUBMED_SCHEMA = deepcopy(DocumentsSchema)
for _column in _PUBMED_SCHEMA.columns.values():
    _column.required = False


_FAILURE_COLUMNS = [
    "schema_context",
    "column",
    "check",
    "check_number",
    "failure_case",
    "index",
]


def _empty_failure_cases() -> pd.DataFrame:
    """Return an empty failure case frame matching pandera layout."""

    return pd.DataFrame(columns=_FAILURE_COLUMNS)


def _build_missing_pmid_failures(series: pd.Series) -> pd.DataFrame:
    """Return failure cases for rows missing ``PubMed.PMID`` values."""

    if series.empty:
        return _empty_failure_cases()

    failures: list[dict[str, object]] = []
    for index, value in series.items():
        failures.append(
            {
                "schema_context": "PubMedDocumentsSchema",
                "column": "PubMed.PMID",
                "check": "not_null",
                "check_number": 0,
                "failure_case": value,
                "index": index,
            }
        )
    return pd.DataFrame(failures, columns=_FAILURE_COLUMNS)


def _validate_documents(df: pd.DataFrame) -> SchemaValidationResult:
    """Validate PubMed document chunks using the relaxed schema."""

    validated = _PUBMED_SCHEMA.validate(df, lazy=True)

    if "PubMed.PMID" not in validated.columns:
        failures = _build_missing_pmid_failures(
            pd.Series(dtype=object, index=validated.index)
        )
        return SchemaValidationResult(
            validated.iloc[0:0], failures, "PubMedDocumentsSchema"
        )

    pmid_series = validated["PubMed.PMID"].astype("string")
    missing_mask = pmid_series.isna() | pmid_series.str.strip().eq("")
    if not missing_mask.any():
        return SchemaValidationResult(
            validated, _empty_failure_cases(), "PubMedDocumentsSchema"
        )

    failure_cases = _build_missing_pmid_failures(pmid_series[missing_mask])
    cleaned = validated.loc[~missing_mask].copy()
    return SchemaValidationResult(cleaned, failure_cases, "PubMedDocumentsSchema")


def _stream_pubmed_batches(
    *,
    pmids: Sequence[str],
    batch_size: int,
    delay: float,
    pubmed_client: PubMedClient,
    limiter: RateLimiter,
    semantic_scholar_cfg: SemanticScholarCfg,
    openalex_cfg: OpenAlexCfg,
    crossref_cfg: CrossRefCfg,
    dump_level: str,
    openalex_limiter: RateLimiter,
    crossref_limiter: RateLimiter,
    api_cfg: ApiCfg,
    retry_cfg: RetryCfg,
) -> Iterator[pd.DataFrame]:
    """Yield combined metadata for ``pmids`` one batch at a time."""

    with session_with_retry(api_cfg, retry_cfg) as session:
        for i in range(0, len(pmids), batch_size):
            batch_pmids = list(pmids[i : i + batch_size])
            limiter.acquire()
            pubmed_list = fetch_pubmed_batch(
                session,
                batch_pmids,
                delay,
                client=pubmed_client,
                retry_cfg=retry_cfg,
            )
            limiter.acquire()
            semsch_list = fetch_semantic_scholar_batch(
                session,
                batch_pmids,
                delay,
                cfg=semantic_scholar_cfg,
                retry_cfg=retry_cfg,
            )
            semsch_map = {s.get("scholar.PMID"): s for s in semsch_list}
            combined_records: list[dict[str, str]] = []
            for pubmed in pubmed_list:
                pmid = pubmed.get("PubMed.PMID", "")
                semsch = semsch_map.get(pmid, {})

                openalex = fetch_openalex(
                    session,
                    pmid,
                    cfg=openalex_cfg,
                    limiter=openalex_limiter,
                    retry_cfg=retry_cfg,
                )
                doi = pubmed.get("PubMed.DOI") or semsch.get("scholar.DOI") or ""
                crossref = fetch_crossref(
                    session,
                    doi,
                    cfg=crossref_cfg,
                    limiter=crossref_limiter,
                    retry_cfg=retry_cfg,
                )

                combined = merge_records(pubmed, semsch, openalex, crossref)
                combined_records.append(combined)

            if combined_records:
                print_results(combined_records, level=dump_level)
                yield pd.DataFrame.from_records(combined_records)


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
        metadata = getattr(args, "_config_metadata", None)
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    if metadata is not None:
        logger.info("config_snapshot", config=getattr(metadata, "snapshot", {}))

    try:
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

    global_rps = cfg.rate.global_rps
    if global_rps is None:
        global_rps = 0
    pubmed_rps = cfg.pubmed.rps if cfg.pubmed.rps is not None else global_rps
    if pubmed_rps is None:
        pubmed_rps = 0
    pubmed_burst = (
        cfg.pubmed.burst if cfg.pubmed.burst is not None else cfg.rate.global_burst
    )
    if pubmed_burst is not None and pubmed_burst <= 0:
        pubmed_burst = None
    limiter = get_limiter("pubmed", pubmed_rps, pubmed_burst)
    delay = 1.0 / pubmed_rps if pubmed_rps > 0 else 0.0

    pmid_df = read_pmids(args.input_csv, cfg=cfg.pubmed)
    pmids = pmid_df["PMID"].tolist()
    openalex_limiter = get_limiter("openalex", cfg.openalex.rps, cfg.openalex.burst)
    crossref_limiter = get_limiter("crossref", cfg.crossref.rps, cfg.crossref.burst)
    pubmed_client = PubMedClient(cfg.pubmed)

    batch_size = cfg.document.pubmed.batch_size
    dump_level = "INFO" if args.keep_verbose_dumps else "DEBUG"

    output_path = (
        Path(args.output_csv)
        if args.output_csv
        else Path(f"output.{Path(args.input_csv).stem}_{date.today():%Y%m%d}.csv")
    )
    failure_path = output_path.with_name(f"{output_path.stem}_failure_cases.csv")

    def fetcher() -> Iterable[pd.DataFrame]:
        return _stream_pubmed_batches(
            pmids=pmids,
            batch_size=batch_size,
            delay=delay,
            pubmed_client=pubmed_client,
            limiter=limiter,
            semantic_scholar_cfg=cfg.semantic_scholar,
            openalex_cfg=cfg.openalex,
            crossref_cfg=cfg.crossref,
            dump_level=dump_level,
            openalex_limiter=openalex_limiter,
            crossref_limiter=crossref_limiter,
            api_cfg=cfg.api,
            retry_cfg=cfg.retry,
        )

    metadata_hooks: list[MetadataHook] = [normalize_documents, add_pipeline_metadata]

    stats_tracker = {"rows_total": 0, "rows_kept": 0}
    validation_failures: list[pd.DataFrame] = []

    def _validator_with_stats(chunk: pd.DataFrame) -> SchemaValidationResult:
        stats_tracker["rows_total"] += len(chunk)
        result = _validate_documents(chunk)
        stats_tracker["rows_kept"] += len(result.data)
        if not result.failure_cases.empty:
            validation_failures.append(result.failure_cases)
            logger.error(
                "validation_failed",
                failures=len(result.failure_cases),
                path=str(failure_path),
            )
        return SchemaValidationResult(
            result.data,
            _empty_failure_cases(),
            "PubMedDocumentsSchema",
        )

    validators: list[Validator] = [cast(Validator, _validator_with_stats)]

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: Sequence[str] | None,
        key_cols: Sequence[str],
    ) -> Path:
        sort_columns = list(key_cols)
        if not sort_columns and col_order is not None:
            sort_columns = list(col_order)
        order = list(col_order) if col_order is not None else sorted(sort_columns)
        return write_csv_chunks_deterministic(
            chunks,
            destination,
            key_cols=sort_columns,
            col_order=order,
            chunksize=cfg.io.csv_chunksize,
            sort_chunksize=cfg.io.csv_chunksize,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
        )

    table_quality = build_table_quality_hook(
        cfg.system.doc_quality,
        table_name=output_path.with_suffix(""),
        destination=output_path.parent,
    )

    command = " ".join(["pubmed_library"] + (list(argv) if argv else []))
    config_snapshot = _serialize_paths(cfg.to_dict())
    inputs = {"input_csv": str(args.input_csv)}
    definition = PipelineDefinition(
        schema=_PUBMED_SCHEMA,
        schema_name="PubMedDocumentsSchema",
        validators=validators,
        metadata_hooks=metadata_hooks,
        writer=writer,
        command=command,
        config_snapshot=config_snapshot,
        inputs=inputs,
        key_columns=["PubMed.PMID"],
        table_quality=table_quality,
    )
    exit_code = run_pipeline(
        definition=definition,
        fetcher=fetcher,
        output_path=output_path,
        failure_path=failure_path,
        cfg=cfg,
        logger=logger,
    )
    if output_path.exists():
        rows_total = stats_tracker["rows_total"]
        rows_kept = stats_tracker["rows_kept"]
        rows_dropped = max(rows_total - rows_kept, 0)
        stats: Stats = {
            "rows_total": rows_total,
            "rows_kept": rows_kept,
            "rows_dropped": rows_dropped,
            "output_sha256": file_sha256(output_path),
        }
        write_meta_yaml(
            csv_path=output_path,
            command=command,
            config_subset=config_snapshot,
            inputs=inputs,
            stats=stats,
            schema="PubMedDocumentsSchema",
        )
    if validation_failures:
        failure_df = pd.concat(validation_failures, ignore_index=True)
        failure_df.to_csv(failure_path, index=False)
        exit_code = 1

    if exit_code == 0:
        logger.info("file_written", path=str(output_path))
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
