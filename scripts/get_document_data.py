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
from collections.abc import Iterable, Mapping, Sequence
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from itertools import chain, islice
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
    positive_int,
)
from library.config import (
    Config,
    CrossRefCfg,
    OpenAlexCfg,
    PubMedCfg,
    SemanticScholarCfg,
    _serialize_paths,
    ensure_dirs,
    print_config,
    session_with_retry,
)
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
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.pipeline_metadata import add_pipeline_metadata
from library.rate_limiter import get_limiter
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from schemas import DocumentsSchema, normalize_documents


def _build_fallback_doi_map(
    frame: pd.DataFrame,
    *,
    pmid_column: str,
    doi_column: str,
) -> dict[str, str]:
    """Return a mapping of PubMed identifiers to DOI overrides."""

    if pmid_column not in frame.columns or doi_column not in frame.columns:
        missing = [
            column
            for column in (pmid_column, doi_column)
            if column not in frame.columns
        ]
        raise ValueError(
            "missing columns in fallback DOI file: "
            f"{', '.join(missing)}"
        )

    mapping: dict[str, str] = {}
    for pmid_value, doi_value in frame[[pmid_column, doi_column]].itertuples(index=False):
        if pd.isna(pmid_value):
            pmid = ""
        else:
            pmid = str(pmid_value).strip()
        if pd.isna(doi_value):
            doi = ""
        else:
            doi = normalise_doi(doi_value)
        if not pmid or not doi:
            continue
        mapping[pmid] = doi
    return mapping


def fetch_pubmed_records(
    pmids: Iterable[str],
    *args: object,
    sleep: float | None = None,
    pubmed_cfg: PubMedCfg | None = None,
    semantic_scholar_cfg: SemanticScholarCfg | None = None,
    openalex_cfg: OpenAlexCfg | None = None,
    crossref_cfg: CrossRefCfg | None = None,
    max_workers: int | None = None,
    batch_size: int | None = None,
    fallback_doi_map: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Retrieve metadata for a sequence of PubMed identifiers.

    Parameters
    ----------
    pmids : Iterable[str]
        Identifiers to query across the PubMed, Semantic Scholar, OpenAlex and
        CrossRef services.
    sleep : float, optional
        Seconds to pause between successive PubMed requests. The same interval
        is reused when the pipeline falls back to the Semantic Scholar
        single-record endpoint, preventing the thread pool from overwhelming
        either API.
    pubmed_cfg : PubMedCfg, optional
        Configuration used when requesting PubMed batches, including base URLs,
        authentication headers and retry behaviour.
    semantic_scholar_cfg : SemanticScholarCfg, optional
        Settings for the Semantic Scholar batch and fallback requests.
    openalex_cfg : OpenAlexCfg, optional
        Connection parameters for OpenAlex lookups.
    crossref_cfg : CrossRefCfg, optional
        Connection parameters for CrossRef lookups driven by the resolved DOI.
    max_workers : int, optional
        Maximum number of concurrent worker threads submitting batches. Higher
        values increase parallelism across both the PubMed batching and the
        Semantic Scholar enrichment stages.
    batch_size : int, optional
        Maximum number of PMIDs per PubMed request. Smaller batches reduce the
        risk of request failures at the cost of more round-trips and additional
        Semantic Scholar fallbacks.
    fallback_doi_map : Mapping[str, str], optional
        Explicit DOI overrides keyed by PMID. When provided, entries supply the
        DOI used to seed downstream CrossRef lookups whenever neither PubMed
        nor Semantic Scholar returns one.

    Returns
    -------
    pandas.DataFrame
        Combined metadata from the different sources.

    Notes
    -----
    For backward compatibility the function also accepts a
    :class:`~library.config.Config` instance as the first positional argument
    after ``pmids``. When supplied, defaults are loaded from
    ``document.pubmed`` and related API sections before calling helper
    functions, and any keyword arguments supplied to
    :func:`fetch_pubmed_records` override those derived settings.

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
        if pubmed_cfg is None:
            pubmed_cfg = cfg.pubmed
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

    settings = cfg or Config()
    rate_cfg = settings.rate
    if pubmed_cfg is None:
        pubmed_cfg = settings.pubmed
    pubmed_limiter = get_limiter("pubmed", rate_cfg.global_rps, rate_cfg.global_burst)
    semantic_limiter = get_limiter(
        "semantic_scholar", rate_cfg.global_rps, rate_cfg.global_burst
    )

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

    def _coerce_batch_argument(*candidates: object) -> list[str]:
        """Return the first iterable batch argument from ``candidates``."""

        for candidate in candidates:
            if isinstance(candidate, Sequence) and not isinstance(
                candidate, str | bytes | bytearray
            ):
                return [str(item) for item in candidate]
        raise TypeError(
            "fetch_pubmed_records() expected a sequence of PMIDs but received"
            f" {tuple(type(c).__name__ for c in candidates)}"
        )

    session_cfg = settings

    def _fetch_batch(
        first: object,
        *rest: object,
        __cfg: Config = session_cfg,
        **__: object,
    ) -> list[dict[str, str]]:
        """Fetch metadata for a batch of PMIDs.

        Each worker opens its own :class:`requests.Session` and retrieves PubMed
        entries for all PMIDs in ``batch`` using a single request. Metadata from
        Semantic Scholar, OpenAlex and CrossRef are then fetched individually
        for each PMID. Exceptions are logged so a failure in one batch does not
        abort the whole process.
        """

        batch_list = _coerce_batch_argument(first, *rest)

        try:

            with session_with_retry(__cfg.api, __cfg.retry) as session:

                pubmed_list = pl.fetch_pubmed_batch(
                    session, batch_list, sleep, cfg=pubmed_cfg
                )

                pmids_in_batch = [p.get("PubMed.PMID", "") for p in pubmed_list]
                semantic_pmids = [pmid for pmid in pmids_in_batch if pmid]

                semsch_map: dict[str, dict[str, str]] = {}
                if semantic_pmids:
                    # Fetch Semantic Scholar data in a single batch
                    semantic_limiter.acquire()
                    semsch_list = ssl.fetch_semantic_scholar_batch(
                        session, semantic_pmids, sleep, cfg=semantic_scholar_cfg
                    )

                    # Create a map for easy lookup
                    semsch_map = {
                        s.get("scholar.PMID"): s for s in semsch_list if s.get("scholar.PMID")
                    }

                    # Fallback to the single-record endpoint when the batch request fails
                    fallback_pmids: list[str] = []
                    seen: set[str] = set()
                    for pmid in semantic_pmids:
                        if pmid in seen:
                            continue
                        seen.add(pmid)
                        record = semsch_map.get(pmid)
                        if record is None or record.get("scholar.Error"):
                            fallback_pmids.append(pmid)
                    for pmid in fallback_pmids:
                        semantic_limiter.acquire()
                        fallback_record = ssl.fetch_semantic_scholar(
                            session, pmid, sleep, cfg=semantic_scholar_cfg
                        )
                        semsch_map[pmid] = fallback_record

                combined_records: list[dict[str, str]] = []
                for pubmed in pubmed_list:
                    pmid = pubmed.get("PubMed.PMID", "")
                    semsch = semsch_map.get(pmid, {}) if pmid else {}

                    # Still fetching these individually for now

                    openalex = ocl.fetch_openalex(
                        session, pmid, openalex_cfg, openalex_limiter
                    )
                    fallback_doi = ""
                    if fallback_doi_map:
                        fallback_doi = fallback_doi_map.get(pmid, "")
                    doi = (
                        pubmed.get("PubMed.DOI")
                        or semsch.get("scholar.DOI")
                        or fallback_doi
                        or ""
                    )
                    crossref = ocl.fetch_crossref(
                        session, doi, crossref_cfg, crossref_limiter
                    )

                    combined = merge_metadata(pubmed, semsch, openalex, crossref)
                    combined_records.append(combined)
                return combined_records
        except requests.RequestException as exc:  # pragma: no cover - network errors
            logger.warning(
                "pubmed_batch_request_failed",
                pmids=batch_list,
                error=str(exc),
            )
            return _failure_records(batch_list, str(exc))
        except Exception as exc:  # pragma: no cover - defensive safety net
            logger.warning(
                "pubmed_batch_unexpected_error",
                pmids=batch_list,
                error=str(exc),
            )
            return _failure_records(batch_list, str(exc))

    iterator = (p for p in pmids if p)
    records: list[dict[str, str]] = []
    tasks: dict[Future[list[dict[str, str]]], tuple[int, list[str]]] = {}
    ordered: dict[int, list[dict[str, str]]] = {}
    processed = 0
    max_in_flight = max(1, max_workers * 2)
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        offset = 0
        pending: set[Future[list[dict[str, str]]]] = set()
        for batch in _chunked(iterator, batch_size):
            if not batch:
                continue
            future = ex.submit(_fetch_batch, batch)
            tasks[future] = (offset, batch)
            pending.add(future)
            offset += len(batch)
            if len(pending) >= max_in_flight:
                done_future = next(as_completed(list(pending)))
                pending.remove(done_future)
                batch_offset, submitted_batch = tasks.pop(done_future)
                ordered[batch_offset] = done_future.result()
                processed += len(submitted_batch)
                logger.info("documents_processed", count=processed)

        while pending:
            done_future = next(as_completed(list(pending)))
            pending.remove(done_future)
            batch_offset, submitted_batch = tasks.pop(done_future)
            ordered[batch_offset] = done_future.result()
            processed += len(submitted_batch)
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


_EXPORT_COLUMNS = [

    "PubMed.PMID",
    "PubMed.DOI",
    "PubMed.ArticleTitle",
    "PubMed.Abstract",
    "PubMed.JournalTitle",
    "PubMed.JournalISOAbbrev",
    "PubMed.Volume",
    "PubMed.Issue",
    "PubMed.StartPage",
    "PubMed.EndPage",
    "PubMed.ISSN",
    "PubMed.PublicationType",
    "PubMed.MeSH_Descriptors",
    "PubMed.MeSH_Qualifiers",
    "PubMed.ChemicalList",
    "PubMed.YearCompleted",
    "PubMed.MonthCompleted",
    "PubMed.DayCompleted",
    "PubMed.YearRevised",
    "PubMed.MonthRevised",
    "PubMed.DayRevised",
    "PubMed.Error",
    "scholar.PMID",
    "scholar.DOI",
    "scholar.PublicationTypes",
    "scholar.Venue",
    "scholar.SemanticScholarId",
    "scholar.ExternalIds",
    "scholar.Error",
    "OpenAlex.PMID",
    "OpenAlex.DOI",
    "OpenAlex.PublicationTypes",
    "OpenAlex.TypeCrossref",
    "OpenAlex.Genre",
    "OpenAlex.Venue",
    "OpenAlex.MeshDescriptors",
    "OpenAlex.MeshQualifiers",
    "OpenAlex.Id",
    "OpenAlex.Error",
    "crossref.DOI",
    "crossref.Type",
    "crossref.Subtype",
    "crossref.Title",
    "crossref.Subtitle",
    "crossref.Subject",
    "crossref.Error",
    "publication_types_normalised",
    "publication_review_score",
    "publication_experimental_score",
    "publication_class",
    "ChEMBL.document_chembl_id",
    "ChEMBL.title",
    "ChEMBL.abstract",
    "ChEMBL.doi",
    "ChEMBL.year",
    "ChEMBL.journal",
    "ChEMBL.journal_abbrev",
    "ChEMBL.volume",
    "ChEMBL.issue",
    "ChEMBL.first_page",
    "ChEMBL.last_page",
    "ChEMBL.pubmed_id",
    "ChEMBL.authors",
    "ChEMBL.source",
]

_EXPORT_COLUMN_RENAMES = {
    "document_chembl_id": "ChEMBL.document_chembl_id",
    "title": "ChEMBL.title",
    "abstract": "ChEMBL.abstract",
    "doi": "ChEMBL.doi",
    "year": "ChEMBL.year",
    "journal": "ChEMBL.journal",
    "journal_abbrev": "ChEMBL.journal_abbrev",
    "volume": "ChEMBL.volume",
    "issue": "ChEMBL.issue",
    "first_page": "ChEMBL.first_page",
    "last_page": "ChEMBL.last_page",
    "pubmed_id": "ChEMBL.pubmed_id",
    "authors": "ChEMBL.authors",
    "source": "ChEMBL.source",
    "publication_type_score_review": "publication_review_score",
    "publication_type_score_experimental": "publication_experimental_score",
}

_EXPORT_COALESCE_SOURCES = {
    "OpenAlex.PMID": ["OpenAlex.PMID", "PubMed.PMID", "scholar.PMID"],
    "OpenAlex.DOI": ["OpenAlex.DOI", "PubMed.DOI", "scholar.DOI", "doi_normalised"],
    "crossref.DOI": ["crossref.DOI", "doi_normalised", "PubMed.DOI", "scholar.DOI"],
}

_EXPORT_SORT_FALLBACK = [
    "ChEMBL.document_chembl_id",
    "PubMed.PMID",
    "scholar.PMID",
    "OpenAlex.PMID",
    "ChEMBL.pubmed_id",
]


def _coalesce_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.Series:
    """Return the first non-empty value across ``columns`` for each row."""

    result = pd.Series("", index=df.index, dtype=object)
    for col in columns:
        if col not in df.columns:
            continue
        values = df[col].fillna("").astype(str)
        mask = result.eq("")
        if mask.any():
            result.loc[mask] = values.loc[mask]
    return result


def _prepare_export_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Rename and project columns to match the export schema."""

    frame = df.copy()
    rename_map = {k: v for k, v in _EXPORT_COLUMN_RENAMES.items() if k in frame.columns}
    if rename_map:
        frame = frame.rename(columns=rename_map)

    for target, sources in _EXPORT_COALESCE_SOURCES.items():
        frame[target] = _coalesce_columns(frame, sources)

    for column in _EXPORT_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""

    return frame[_EXPORT_COLUMNS]


def _finalise_export(
    df: pd.DataFrame,
    output: Path,
    cfg: Config,
    *,
    input_csv: Path,
    key_columns: Sequence[str] | None = None,
) -> int:
    """Validate ``df`` and write CSV/metadata artefacts."""

    df = add_pipeline_metadata(df)
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
                "missing_optional_columns",
                columns=sorted(missing_optional),
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
                "document_validation_failed",
                failure_count=len(exc.failure_cases),
                failure_path=str(failure_path),
                error=str(exc),
            )
            validated = getattr(exc, "validated_data", ordered)
            exit_code = 1
    else:
        logger.warning(
            "validation_skipped_missing_required",
            columns=sorted(missing_required),
        )

    rows_kept = len(validated)
    rows_dropped = rows_total - rows_kept

    export_df = build_dataframe(
        validated, columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False
    )
    export_df = dataframe_to_strings(export_df, skip=_NUMERIC_EXPORT_COLUMNS)
    export_df = _prepare_export_frame(export_df)

    key_cols: list[str] = []
    if key_columns:
        for column in key_columns:
            mapped = _EXPORT_COLUMN_RENAMES.get(column, column)
            if mapped in export_df.columns and mapped not in key_cols:
                key_cols.append(mapped)
    if not key_cols:
        for candidate in _EXPORT_SORT_FALLBACK:
            if candidate in export_df.columns:
                key_cols = [candidate]
                break
    if not key_cols:
        key_cols = [_EXPORT_COLUMNS[0]]

    col_order = list(_EXPORT_COLUMNS)
    try:
        csv_path = io.write_csv(
            export_df,
            output,
            cfg=cfg,
            key_cols=key_cols,
            col_order=col_order,
        )
    except OSError as exc:
        logger.error("csv_write_failed", error=str(exc), path=str(output))
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

    quality_path = csv_path.with_suffix(".quality.json")
    try:
        report = build_quality_report(validated)
        save_quality_report(report, quality_path)
    except (OSError, TypeError, ValueError) as exc:
        logger.error(
            "quality_report_write_failed",
            error=str(exc),
            path=str(quality_path),
        )
        return 1

    try:
        analyze_table_quality(validated, table_name=str(csv_path.with_suffix("")))
    except ValueError as exc:
        logger.error("quality_report_generation_failed", error=str(exc))
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
        logger.error(
            "invalid_limit",
            section="document.pubmed",
            limit=limit,
        )
        return 1
    try:
        pmids_iter = io.read_ids(
            args.input_csv, column=cfg.document.pubmed.column, cfg=cfg.io
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "input_read_failed",
            error=str(exc),
            path=str(args.input_csv),
        )
        return 1
    pmids: Iterable[str] = pmids_iter
    if limit is not None:
        limited_pmids = list(islice(pmids_iter, limit))
        pmids = limited_pmids
        logger.info("process_limit", limit=len(limited_pmids))

    fallback_doi_map: Mapping[str, str] | None = None
    fallback_csv = getattr(args, "fallback_doi_csv", None)
    if fallback_csv:
        try:
            fallback_frame = io.read_csv(fallback_csv, cfg=cfg.io)
        except (FileNotFoundError, ValueError) as exc:
            logger.error(
                "fallback_doi_read_failed",
                error=str(exc),
                path=str(fallback_csv),
            )
            return 1
        try:
            fallback_map = _build_fallback_doi_map(
                fallback_frame,
                pmid_column=getattr(args, "fallback_doi_pmid_column", "PMID"),
                doi_column=getattr(args, "fallback_doi_value_column", "DOI"),
            )
        except ValueError as exc:
            logger.error(
                "fallback_doi_invalid",
                error=str(exc),
                path=str(fallback_csv),
            )
            return 1
        fallback_doi_map = fallback_map or None

    try:
        df = fetch_pubmed_records(
            pmids,
            cfg,
            sleep=cfg.document.pubmed.sleep,
            pubmed_cfg=cfg.pubmed,
            semantic_scholar_cfg=cfg.semantic_scholar,
            openalex_cfg=cfg.openalex,
            crossref_cfg=cfg.crossref,
            max_workers=cfg.document.pubmed.workers,
            batch_size=cfg.document.pubmed.batch_size,
            fallback_doi_map=fallback_doi_map,
        )
        output = Path(args.output_csv or io.default_output_path(args.input_csv, cfg.io))
        df = normalize_documents(df)
        exit_code = _finalise_export(
            df,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("pubmed_pipeline_failed", error=str(exc))
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
        logger.error("invalid_limit", section="document.chembl", limit=limit)
        return 1

    # Configure session for ChEMBL requests
    with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
        try:
            ids_iter = io.read_ids(
                args.input_csv, column=cfg.document.chembl.column, cfg=cfg.io
            )
        except (FileNotFoundError, ValueError) as exc:
            logger.error(
                "input_read_failed",
                error=str(exc),
                path=str(args.input_csv),
            )
            return 1

        ids: Iterable[str] = ids_iter
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
            logger.error(
                "chembl_documents_fetch_failed",
                error=str(exc),
                chunk_size=cfg.document.chembl.chunk_size,
            )
            return 1
        if "doi" in df.columns:
            df["doi"] = df["doi"].map(normalise_doi)
        output = Path(args.output_csv or io.default_output_path(args.input_csv, cfg.io))
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
        logger.error("invalid_limit", section="document.all", limit=limit)
        return 1

    # Prepare shared session before performing any API calls
    try:
        ids_iter = io.read_ids(
            args.input_csv, column=cfg.document.all.column, cfg=cfg.io
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "input_read_failed",
            error=str(exc),
            path=str(args.input_csv),
        )
        return 1

    ids_source: Iterable[str] = ids_iter
    if limit is not None:
        limited_ids = list(islice(ids_iter, limit))
        ids_source = limited_ids
        logger.info("process_limit", limit=len(limited_ids))

    iterator = iter(ids_source)
    sample_size = cfg.document.all.chunk_size
    sample_ids = list(islice(iterator, sample_size))
    ids_for_fetch = chain(sample_ids, iterator)
    try:
        with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
            doc_df = cl.get_documents(
                ids_for_fetch,
                cfg=cfg.api,
                client=client,
                chunk_size=cfg.document.all.chunk_size,
                timeout=cfg.document.all.timeout,
            )
    except (requests.RequestException, ValueError) as exc:
        logger.error(
            "chembl_documents_fetch_failed",
            ids=sample_ids,
            error=str(exc),
            chunk_size=cfg.document.all.chunk_size,
        )
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
    doi_fallback_map: dict[str, str] = {}
    if "doi" in doc_df.columns:
        doi_series = doc_df["doi"].astype("string")
        mask = pubmed_ids.notna() & doi_series.notna() & (doi_series != "")
        if mask.any():
            masked_pmids = pubmed_ids[mask].tolist()
            masked_dois = doi_series[mask].tolist()
            if len(masked_pmids) != len(masked_dois):
                raise ValueError("mismatched DOI fallback map lengths")
            doi_fallback_map = {
                str(pmid): str(doi)
                for pmid, doi in zip(masked_pmids, masked_dois, strict=True)
            }
    pmids = pubmed_ids.dropna().astype(str).tolist()
    pub_df = fetch_pubmed_records(
        pmids,
        cfg.document.all.sleep,
        cfg.semantic_scholar,
        cfg.openalex,
        cfg.crossref,
        cfg.document.all.workers,
        cfg.document.all.batch_size,
        pubmed_cfg=cfg.pubmed,
        fallback_doi_map=doi_fallback_map or None,
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
        type=Path,
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
        logger.error("config_override_error", error=str(exc))
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("directory_setup_failed", error=str(exc))
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
