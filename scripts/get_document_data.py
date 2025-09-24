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

    python scripts/get_document_data.py pubmed --config config.yaml --input pmids.csv --output output.csv

The input file must contain a ``PMID`` column.

"""

from __future__ import annotations

import sys

# ruff: noqa: E402
from pathlib import Path

if __package__ is None:  # running as a script
    sys.path.append(str(Path(__file__).resolve().parents[1]))

import argparse
from collections.abc import Iterable, Sequence
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from itertools import islice
from typing import cast

import pandas as pd
import requests
from pandera.errors import SchemaErrors

from library import chembl_library as cl
from library import document_postprocessing as dp
from library import io
from library import openalex_crossref_library as ocl
from library import pubmed_library as pl
from library import semantic_scholar_library as ssl
from library.chembl_client import ChemblClient, _chunked
from library.cli import (
    LoggerConfig,
    apply_config_overrides,
    build_root_parser,
    configure_logger,
)
from library.config import (
    ApiCfg,
    Config,
    CrossRefCfg,
    OpenAlexCfg,

    PubMedCfg,
    RetryCfg,
    SemanticScholarCfg,
    _serialize_paths,
    ensure_dirs,
    print_config,
    session_with_retry,
)
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.rate_limiter import get_limiter
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from library.document_pipeline import (
    DOCUMENT_SCHEMA_COLUMNS,
    build_dataframe,
    build_quality_report,
    dataframe_to_strings,
    merge_metadata,
    merge_with_chembl,
    normalise_doi,
    save_quality_report,
)
from schemas import DocumentsSchema, normalize_documents


def fetch_pubmed_records(
    pmids: Iterable[str],

    *args: object,
    sleep: float | None = None,
    semantic_scholar_cfg: SemanticScholarCfg | None = None,
    openalex_cfg: OpenAlexCfg | None = None,
    crossref_cfg: CrossRefCfg | None = None,
    max_workers: int | None = None,
    batch_size: int | None = None,

) -> pd.DataFrame:
    """Retrieve metadata for a sequence of PubMed identifiers.

    Parameters
    ----------
    pmids:
        Identifiers to query.
    sleep:

        Seconds to pause between PubMed and Semantic Scholar requests.
    semantic_scholar_cfg:
        Configuration for Semantic Scholar API access.
    openalex_cfg:
        Configuration for OpenAlex API access.
    crossref_cfg:
        Configuration for CrossRef API access.
    api_cfg:
        Shared API configuration providing retry defaults and headers.
    retry_cfg:
        Retry behaviour applied to outbound HTTP sessions.
    pubmed_cfg:
        Settings for the PubMed API, including timeouts and retry counts.
    semantic_cfg:
        Semantic Scholar configuration used for batch enrichment.
    max_workers:
        Maximum number of concurrent threads.
    batch_size:
        Maximum number of PMIDs per PubMed request.

    Returns
    -------
    pandas.DataFrame
        Combined metadata from the different sources.

    Notes
    -----
    For backward compatibility the function also accepts a
    :class:`~library.config.Config` instance as the first positional argument
    after ``pmids``. When supplied, connection parameters and batching options
    are derived from ``config.document.pubmed`` unless overridden explicitly
    via keyword arguments.

    """

    cfg: Config | None = None
    positional = list(args)
    if positional:
        candidate = positional[0]
        if isinstance(candidate, Config):
            cfg = candidate
            positional = positional[1:]
        elif sleep is None:
            if len(positional) < 4:
                raise TypeError(
                    "fetch_pubmed_records() missing required positional arguments: "
                    "'sleep', 'semantic_scholar_cfg', 'openalex_cfg', 'crossref_cfg'"
                )
            sleep = cast(float, positional[0])
            semantic_scholar_cfg = cast(SemanticScholarCfg, positional[1])
            openalex_cfg = cast(OpenAlexCfg, positional[2])
            crossref_cfg = cast(CrossRefCfg, positional[3])
            if len(positional) > 4:
                max_workers = cast(int, positional[4])
            if len(positional) > 5:
                batch_size = cast(int, positional[5])
            positional = []
        else:
            raise TypeError(
                "fetch_pubmed_records() received multiple values for 'sleep'"
            )
    if positional:
        raise TypeError("fetch_pubmed_records() got unexpected positional arguments")

    if cfg is not None:
        if sleep is None:
            sleep = cfg.document.pubmed.sleep
        if semantic_scholar_cfg is None:
            semantic_scholar_cfg = cfg.semantic_scholar
        if openalex_cfg is None:
            openalex_cfg = cfg.openalex
        if crossref_cfg is None:
            crossref_cfg = cfg.crossref
        if max_workers is None:
            max_workers = cfg.document.pubmed.workers
        if batch_size is None:
            batch_size = cfg.document.pubmed.batch_size

    if (
        sleep is None
        or semantic_scholar_cfg is None
        or openalex_cfg is None
        or crossref_cfg is None
    ):
        raise TypeError(
            "fetch_pubmed_records() missing required configuration. "
            "Provide either a Config instance as the second positional argument "
            "or explicit keyword arguments."
        )

    if max_workers is None:
        max_workers = 1
    if batch_size is None:
        batch_size = 100
    openalex_limiter = get_limiter("openalex", openalex_cfg.rps, openalex_cfg.burst)
    crossref_limiter = get_limiter("crossref", crossref_cfg.rps, crossref_cfg.burst)

    def _failure_records(batch: Sequence[str], message: str) -> list[dict[str, str]]:
        """Return placeholder rows describing a failure for ``batch`` PMIDs."""

        records: list[dict[str, str]] = []
        for pmid in batch:
            pubmed = {"PubMed.PMID": pmid, "PubMed.Error": message}
            scholar = {"scholar.PMID": pmid, "scholar.Error": message}
            openalex = {"OpenAlex.Error": message}
            crossref = {"crossref.Error": message}
            records.append(merge_metadata(pubmed, scholar, openalex, crossref))
        return records

    def _fetch_batch(batch: list[str]) -> list[dict[str, str]]:

        """Fetch metadata for a batch of PMIDs.

        Each worker opens its own :class:`requests.Session` and retrieves PubMed
        entries for all PMIDs in ``batch`` using a single request. Metadata from
        Semantic Scholar, OpenAlex and CrossRef are then fetched individually
        for each PMID. Exceptions are logged so a failure in one batch does not
        abort the whole process.
        """
        try:
            with session_with_retry(api_cfg, retry_cfg) as session:
                pubmed_limiter.acquire()
                pubmed_list = pl.fetch_pubmed_batch(
                    session, batch, sleep, cfg=pubmed_cfg
                )
                pmids_in_batch = [p.get("PubMed.PMID", "") for p in pubmed_list]

                # Fetch Semantic Scholar data in a single batch
                semantic_limiter.acquire()
                semsch_list = ssl.fetch_semantic_scholar_batch(
                    session, pmids_in_batch, sleep, cfg=semantic_scholar_cfg

                )

                # Create a map for easy lookup
                semsch_map = {s.get("scholar.PMID"): s for s in semsch_list}

                combined_records: list[dict[str, str]] = []
                for pubmed in pubmed_list:
                    pmid = pubmed.get("PubMed.PMID", "")
                    semsch = semsch_map.get(pmid, {})

                    # Still fetching these individually for now

                    openalex = ocl.fetch_openalex(
                        session, pmid, openalex_cfg, openalex_limiter
                    )
                    doi = pubmed.get("PubMed.DOI") or semsch.get("scholar.DOI") or ""
                    crossref = ocl.fetch_crossref(
                        session, doi, crossref_cfg, crossref_limiter
                    )

                    combined = merge_metadata(pubmed, semsch, openalex, crossref)
                    combined_records.append(combined)
                return combined_records
        except requests.RequestException as exc:  # pragma: no cover - network errors
            logger.warning("failed to fetch PMIDs %s: %s", batch, exc)
            return _failure_records(batch, str(exc))

        except Exception as exc:  # pragma: no cover - defensive safety net

            logger.warning("unexpected error for PMIDs %s: %s", batch, exc)
            return _failure_records(batch, str(exc))

    iterator = (p for p in pmids if p)
    records: list[dict[str, str]] = []
    tasks: dict[Future[list[dict[str, str]]], tuple[int, list[str]]] = {}
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        offset = 0
        for batch in _chunked(iterator, batch_size):
            if not batch:
                continue
            future = ex.submit(_fetch_batch, offset, batch)
            tasks[future] = (offset, batch)
            offset += len(batch)

        processed = 0
        ordered: dict[int, list[dict[str, str]]] = {}
        for future in as_completed(tasks):
            offset, batch = tasks[future]
            result = future.result()
            ordered[offset] = result
            processed += len(batch)
            logger.info("documents_processed", count=processed)

    for offset in sorted(ordered):
        records.extend(ordered[offset])
    if not records:
        return build_dataframe([], columns=DOCUMENT_SCHEMA_COLUMNS)
    return build_dataframe(records, columns=DOCUMENT_SCHEMA_COLUMNS)


_NUMERIC_EXPORT_COLUMNS = {
    "publication_type_score_review",
    "publication_type_score_experimental",
    "publication_type_score_unknown",
}


def _finalise_export(
    df: pd.DataFrame,
    output: Path,
    cfg: Config,
    *,
    input_csv: Path,
    key_columns: Sequence[str] | None = None,
) -> int:
    """Validate ``df`` and write CSV/metadata artefacts."""

    ordered = build_dataframe(df, columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False)
    rows_total = len(ordered)
    exit_code = 0
    required_cols = {
        name for name, col in DocumentsSchema.columns.items() if col.required
    }
    optional_cols = set(DocumentsSchema.columns) - required_cols
    missing_required = required_cols - set(ordered.columns)
    missing_optional = optional_cols - set(ordered.columns)

    validated = ordered
    if not missing_required:
        if missing_optional:
            logger.warning(
                "DataFrame is missing optional columns: %s", missing_optional
            )
        try:
            validated = DocumentsSchema.validate(ordered, lazy=True)
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
            validated = getattr(exc, "validated_data", ordered)
            exit_code = 1
    else:
        logger.warning(
            "Skipping validation due to missing required columns: %s",
            missing_required,
        )

    rows_kept = len(validated)
    rows_dropped = rows_total - rows_kept

    export_df = build_dataframe(
        validated, columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False
    )
    export_df = dataframe_to_strings(export_df, skip=_NUMERIC_EXPORT_COLUMNS)
    if key_columns:
        key_cols = [c for c in key_columns if c in export_df.columns]
    else:
        key_cols = []
    col_order = [
        c for c in DOCUMENT_SCHEMA_COLUMNS if c in export_df.columns
    ] + sorted(c for c in export_df.columns if c not in DOCUMENT_SCHEMA_COLUMNS)
    try:
        csv_path = io.write_csv(
            export_df,
            output,
            cfg=cfg,
            key_cols=key_cols or None,
            col_order=col_order,
        )
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1
    logger.info("write_done", rows=rows_kept, path=str(csv_path))

    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_kept,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(csv_path),
    }
    write_meta_yaml(
        csv_path=csv_path,
        command=" ".join(sys.argv),
        config_subset=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(input_csv)},
        stats=stats,
        schema="DocumentsSchema",
    )

    try:
        report = build_quality_report(validated)
        save_quality_report(report, csv_path.with_suffix(".quality.json"))
    except (OSError, TypeError, ValueError) as exc:
        logger.error("failed to write quality report: %s", exc)
        return 1

    try:
        analyze_table_quality(validated, table_name=str(csv_path.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return exit_code


def run_pubmed(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``pubmed`` sub-command.

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
    limit = cfg.document.pubmed.limit
    if limit is not None and limit < 0:
        logger.error("document.pubmed.limit must be non-negative")
        return 1
    try:
        pmids_iter = io.read_ids(
            args.input_csv, column=cfg.document.pubmed.column, cfg=cfg.io
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1
    pmids: Iterable[str] = pmids_iter
    if limit is not None:
        limited_pmids = list(islice(pmids_iter, limit))
        pmids = limited_pmids
        logger.info("process_limit", limit=len(limited_pmids))

    try:
        df = fetch_pubmed_records(
            pmids,

            cfg.document.pubmed.sleep,
            cfg.semantic_scholar,
            cfg.openalex,
            cfg.crossref,
            cfg.document.pubmed.workers,
            cfg.document.pubmed.batch_size,
        )
        output = Path(
            args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        )

        df = normalize_documents(df)
        exit_code = _finalise_export(
            df,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1
    return exit_code


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``chembl`` sub-command.

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
    limit = cfg.document.chembl.limit
    if limit is not None and limit < 0:
        logger.error("document.chembl.limit must be non-negative")
        return 1

    # Configure session for ChEMBL requests
    with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
        try:
            ids_iter = io.read_ids(
                args.input_csv, column=cfg.document.chembl.column, cfg=cfg.io
            )
        except (FileNotFoundError, ValueError) as exc:
            logger.error("%s", exc)
            return 1

        ids = ids_iter
        if limit is not None:
            limited_ids = list(islice(ids_iter, limit))
            ids = limited_ids
            logger.info("process_limit", limit=len(limited_ids))

        try:
            df = cl.get_documents(
                ids,
                cfg=cfg.api,
                client=client,
                chunk_size=cfg.document.chembl.chunk_size,
                timeout=cfg.document.chembl.timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error("failed to retrieve documents: %s", exc)
            return 1
        if "doi" in df.columns:
            df["doi"] = df["doi"].map(normalise_doi)
        output = Path(
            args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        )
        df = normalize_documents(df)
        exit_code = _finalise_export(
            df,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
        )
        return exit_code


def run_all(cfg: Config, args: argparse.Namespace) -> int:
    """Run ChEMBL and PubMed pipelines and merge their outputs.

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
    limit = cfg.document.all.limit
    if limit is not None and limit < 0:
        logger.error("document.all.limit must be non-negative")
        return 1

    # Prepare shared session before performing any API calls
    try:
        ids_iter = io.read_ids(
            args.input_csv, column=cfg.document.all.column, cfg=cfg.io
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    ids_source: Iterable[str] = ids_iter
    if limit is not None:
        limited_ids = list(islice(ids_iter, limit))
        ids_source = limited_ids
        logger.info("process_limit", limit=len(limited_ids))

    chunk_ids: list[str] = []
    try:
        with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
            chunk_ids = list(ids_source)  # capture IDs for logging on failure
            doc_df = cl.get_documents(
                chunk_ids,
                cfg=cfg.api,
                client=client,
                chunk_size=cfg.document.all.chunk_size,
                timeout=cfg.document.all.timeout,
            )
    except (requests.RequestException, ValueError) as exc:
        logger.error("failed to retrieve documents for %s: %s", chunk_ids, exc)
        return 1
    output = Path(args.output_csv or io.default_output_path(args.input_csv, cfg.io))
    if "doi" in doc_df.columns:
        doc_df["doi"] = doc_df["doi"].map(normalise_doi)
    if doc_df.empty or "pubmed_id" not in doc_df:
        processed = dp.postprocess_documents(doc_df)
        extra_cols = [c for c in doc_df.columns if c not in processed.columns]
        if extra_cols:
            processed = processed.merge(
                doc_df[["document_chembl_id"] + extra_cols],
                on="document_chembl_id",
                how="left",
            )
        processed = normalize_documents(processed)
        exit_code = _finalise_export(
            processed,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
        )
        return exit_code

    # Normalise PubMed identifiers to strings to avoid dtype mismatches
    pubmed_ids = pd.to_numeric(doc_df["pubmed_id"], errors="coerce").astype("Int64")
    pmids = pubmed_ids.dropna().astype(str).tolist()
    pub_df = fetch_pubmed_records(
        pmids,

        cfg.document.all.sleep,
        cfg.semantic_scholar,
        cfg.openalex,
        cfg.crossref,
        cfg.document.all.workers,
        cfg.document.all.batch_size,
    )
    doc_df["pubmed_id"] = pubmed_ids.astype("Int64").astype("string").fillna("")
    merged = merge_with_chembl(doc_df, pub_df)
    processed = dp.postprocess_documents(merged)
    extra_cols = [c for c in merged.columns if c not in processed.columns]
    if extra_cols:
        processed = processed.merge(
            merged[["document_chembl_id"] + extra_cols],
            on="document_chembl_id",
            how="left",
        )
    processed = normalize_documents(processed)
    exit_code = _finalise_export(
        processed,
        output,
        cfg,
        input_csv=Path(args.input_csv),
        key_columns=["document_chembl_id"],
    )
    return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the argument parser for document utilities.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        Parser populated with all sub-commands and logging configuration.

    """
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
        type=int,
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
        "--chunk-size", type=int, default=5, help="Maximum number of IDs per request"
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
        "--chunk-size", type=int, default=5, help="Maximum IDs per request"
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
        type=int,
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
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
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
        cfg: Config = apply_config_overrides(
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
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error("%s", exc)
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("failed to set up directories: %s", exc)
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    exit_code = cast(int, args.func(cfg, args))
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
