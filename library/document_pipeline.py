"""Utilities for constructing the document metadata pipeline.

The helpers defined here provide a central place for declaring column
ordering, normalising metadata fields across sources and producing summary
reports.  They mirror the behaviour of the production pipeline implemented in
the downstream repository and allow :mod:`scripts.get_document_data` to share a
consistent implementation across its ``pubmed``, ``chembl`` and ``all``
sub-commands.
"""

from __future__ import annotations

import json
import argparse
import sys
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from contextlib import ExitStack
from itertools import chain, islice, tee
from pathlib import Path
from typing import Any, cast

import pandas as pd
import requests

from library import chembl_library as cl
from library import document_postprocessing as dp
from library import io
from library import openalex_crossref_library as ocl
from library import pubmed_library as pl
from library import semantic_scholar_library as ssl
from library.clients import ChemblClient, _chunked
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
from library.csv_utils import write_csv_chunks_deterministic
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.pipeline_metadata import add_pipeline_metadata
from library.rate_limiter import RateLimiter, get_limiter
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality

from .chembl_document import DOCUMENT_COLUMNS as _CHEMBL_COLUMNS
from .document_type_classifier import compute_scores, decide_label
from .document_type_terms import parse_terms

# ---------------------------------------------------------------------------
# Column declarations

# Keep local copies so callers cannot mutate the original definitions imported
# from :mod:`library.chembl_document`.
CH_EMBL_COLUMNS: list[str] = list(_CHEMBL_COLUMNS)

PUBMED_COLUMNS: list[str] = [
    "PubMed.PMID",
    "PubMed.DOI",
    "PubMed.ArticleTitle",
    "PubMed.Abstract",
    "PubMed.JournalTitle",
    "PubMed.Volume",
    "PubMed.Issue",
    "PubMed.StartPage",
    "PubMed.EndPage",
    "PubMed.PublicationType",
    "PubMed.MeSH_Descriptors",
    "PubMed.MeSH_Qualifiers",
    "PubMed.ChemicalList",
    "PubMed.DayRevised",
    "PubMed.MonthRevised",
    "PubMed.YearRevised",
    "PubMed.YearCompleted",
    "PubMed.MonthCompleted",
    "PubMed.DayCompleted",
    "PubMed.Error",
    "PubMed.ISSN",
]

SEMANTIC_SCHOLAR_COLUMNS: list[str] = [
    "scholar.PMID",
    "scholar.Venue",
    "scholar.PublicationTypes",
    "scholar.SemanticScholarId",
    "scholar.ExternalIds",
    "scholar.DOI",
    "scholar.Error",
]

OPENALEX_COLUMNS: list[str] = [
    "OpenAlex.PublicationTypes",
    "OpenAlex.TypeCrossref",
    "OpenAlex.Genre",
    "OpenAlex.Id",
    "OpenAlex.Venue",
    "OpenAlex.MeshDescriptors",
    "OpenAlex.MeshQualifiers",
    "OpenAlex.Error",
]

CROSSREF_COLUMNS: list[str] = [
    "crossref.Type",
    "crossref.Subtype",
    "crossref.Title",
    "crossref.Subtitle",
    "crossref.Subject",
    "crossref.Error",
]

DERIVED_COLUMNS: list[str] = [
    "doi_normalised",
    "publication_types_normalised",
    "publication_type_score_review",
    "publication_type_score_experimental",
    "publication_type_score_unknown",
    "publication_class",
]

# The combined list intentionally mirrors the downstream document pipeline so
# CSV exports preserve a deterministic order.
DOCUMENT_SCHEMA_COLUMNS: list[str] = (
    CH_EMBL_COLUMNS
    + DERIVED_COLUMNS
    + PUBMED_COLUMNS
    + SEMANTIC_SCHOLAR_COLUMNS
    + OPENALEX_COLUMNS
    + CROSSREF_COLUMNS
    + [
        "date_code",
        "Index",
        "PubMed.is_review",
        "scholar.is_review",
        "OpenAlex.is_review",
        "pipeline_version",
        "timestamp_utc",
    ]
)

# Remove accidental duplicates while preserving declaration order. Downstream
# code assumes a one-to-one mapping between column name and position.
DOCUMENT_SCHEMA_COLUMNS = list(dict.fromkeys(DOCUMENT_SCHEMA_COLUMNS))


# ---------------------------------------------------------------------------
# Normalisation helpers


def normalise_doi(value: Any) -> str:
    """Return a canonical DOI string or ``""`` if ``value`` is falsy."""

    if value is None:
        return ""
    doi = str(value).strip()
    if not doi:
        return ""
    doi = doi.lower()
    for prefix in ("doi:", "https://doi.org/", "http://doi.org/", "doi.org/"):
        if doi.startswith(prefix):
            doi = doi[len(prefix) :].strip()
    return doi


def _collect_terms(record: Mapping[str, Any]) -> list[str]:
    """Extract publication type terms from ``record``."""

    terms: list[str] = []
    for key in (
        "PubMed.PublicationType",
        "scholar.PublicationTypes",
        "OpenAlex.PublicationTypes",
        "OpenAlex.TypeCrossref",
        "OpenAlex.Genre",
        "crossref.Type",
        "crossref.Subtype",
    ):
        value = record.get(key)
        terms.extend(parse_terms(value))
    seen: set[str] = set()
    unique_terms: list[str] = []
    for term in terms:
        if term not in seen:
            seen.add(term)
            unique_terms.append(term)
    return sorted(unique_terms)


def merge_metadata(*records: Mapping[str, Any]) -> dict[str, Any]:
    """Merge metadata and compute derived document fields."""

    merged: dict[str, Any] = {}
    for rec in records:
        merged.update(rec)

    # Normalise DOI fields and prefer PubMed over Semantic Scholar when
    # deriving the canonical DOI value.
    pubmed_doi = normalise_doi(merged.get("PubMed.DOI"))
    scholar_doi = normalise_doi(merged.get("scholar.DOI"))
    merged["PubMed.DOI"] = pubmed_doi
    merged["scholar.DOI"] = scholar_doi
    merged["doi"] = normalise_doi(merged.get("doi")) or pubmed_doi or scholar_doi
    merged["doi_normalised"] = merged["doi"]

    # Publication type normalisation combines terms from all sources to ensure
    # deterministic ordering.
    all_terms = _collect_terms(merged)
    merged["publication_types_normalised"] = "; ".join(all_terms)

    # Score and classify the publication based on the available terms.
    pubmed_terms = parse_terms(merged.get("PubMed.PublicationType"))
    scholar_terms = parse_terms(merged.get("scholar.PublicationTypes"))
    openalex_terms = parse_terms(merged.get("OpenAlex.PublicationTypes"))
    openalex_terms += parse_terms(merged.get("OpenAlex.TypeCrossref"))
    openalex_terms += parse_terms(merged.get("OpenAlex.Genre"))
    openalex_terms += parse_terms(merged.get("crossref.Type"))
    openalex_terms += parse_terms(merged.get("crossref.Subtype"))

    scores = compute_scores(pubmed_terms, scholar_terms, openalex_terms)
    merged["publication_type_score_review"] = scores["review"]
    merged["publication_type_score_experimental"] = scores["experimental"]
    merged["publication_type_score_unknown"] = scores["unknown"]
    merged["publication_class"] = decide_label(scores)

    return merged


# ---------------------------------------------------------------------------
# DataFrame utilities


def build_dataframe(
    data: Sequence[Mapping[str, Any]] | pd.DataFrame,
    *,
    columns: Sequence[str] = DOCUMENT_SCHEMA_COLUMNS,
    fill_missing: bool = True,
) -> pd.DataFrame:
    """Return a :class:`~pandas.DataFrame` with deterministic column order."""

    unique_columns = list(dict.fromkeys(columns))

    if isinstance(data, pd.DataFrame):
        df = data.copy()
    else:
        if not data:
            return pd.DataFrame(columns=unique_columns)
        df = pd.DataFrame.from_records(data)

    if fill_missing:
        for col in unique_columns:
            if col not in df.columns:
                df[col] = ""

    extra = sorted(c for c in df.columns if c not in unique_columns)
    head = [c for c in unique_columns if c in df.columns]
    ordered = head + extra if not fill_missing else unique_columns + extra
    return df[ordered]


def merge_with_chembl(
    chembl_df: pd.DataFrame, metadata_df: pd.DataFrame
) -> pd.DataFrame:
    """Merge PubMed style metadata into ``chembl_df``."""

    if metadata_df.empty or "PubMed.PMID" not in metadata_df.columns:
        return chembl_df.copy()

    left = chembl_df.copy()
    right = metadata_df.copy()

    drop_cols = [col for col in CH_EMBL_COLUMNS if col in right.columns]
    if drop_cols:
        right = right.drop(columns=drop_cols)

    if "pubmed_id" in left.columns:
        left["pubmed_id"] = (
            pd.to_numeric(left["pubmed_id"], errors="coerce")
            .astype("Int64")
            .astype("string")
            .fillna("")
        )

    right["PubMed.PMID"] = (
        pd.to_numeric(right["PubMed.PMID"], errors="coerce")
        .astype("Int64")
        .astype("string")
        .fillna("")
    )

    merged = left.merge(
        right,
        how="left",
        left_on="pubmed_id" if "pubmed_id" in left.columns else "PubMed.PMID",
        right_on="PubMed.PMID",
    )
    return merged


def dataframe_to_strings(
    df: pd.DataFrame, *, skip: Iterable[str] | None = None
) -> pd.DataFrame:
    """Convert non-numeric columns of ``df`` to strings."""

    skip_set = set(skip or [])
    result = df.copy()
    for col in result.columns:
        if col in skip_set:
            continue
        series = result[col]
        if pd.api.types.is_numeric_dtype(series.dtype):
            continue
        result[col] = series.astype(str)
    return result


# ---------------------------------------------------------------------------
# Quality reporting


def build_quality_report(df: pd.DataFrame) -> dict[str, Any]:
    """Return a JSON serialisable quality summary for ``df``."""

    total = int(len(df))

    def _coverage(series: pd.Series) -> float:
        if series.empty:
            return 0.0
        mask = series.astype(str).str.strip().astype(bool)
        return float(mask.mean())

    doi_series = (
        df["doi"]
        if "doi" in df.columns
        else df.get("doi_normalised", pd.Series(dtype=str))
    )
    class_counts = (
        df.get("publication_class", pd.Series(dtype=str))
        .fillna("unknown")
        .astype(str)
        .str.strip()
        .replace("", "unknown")
        .value_counts()
        .to_dict()
    )

    error_counts = {
        "pubmed": int(
            df.get("PubMed.Error", pd.Series(dtype=str))
            .astype(str)
            .str.strip()
            .astype(bool)
            .sum()
        ),
        "semantic_scholar": int(
            df.get("scholar.Error", pd.Series(dtype=str))
            .astype(str)
            .str.strip()
            .astype(bool)
            .sum()
        ),
        "openalex": int(
            df.get("OpenAlex.Error", pd.Series(dtype=str))
            .astype(str)
            .str.strip()
            .astype(bool)
            .sum()
        ),
        "crossref": int(
            df.get("crossref.Error", pd.Series(dtype=str))
            .astype(str)
            .str.strip()
            .astype(bool)
            .sum()
        ),
    }

    return {
        "rows_total": total,
        "doi_coverage": _coverage(doi_series),
        "publication_class_counts": class_counts,
        "error_counts": error_counts,
    }


def save_quality_report(report: Mapping[str, Any], path: Path) -> Path:
    """Persist ``report`` as JSON file and return ``path``."""

    path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Export helpers


NUMERIC_EXPORT_COLUMNS: set[str] = {
    "publication_type_score_review",
    "publication_type_score_experimental",
    "publication_type_score_unknown",
}


DOCUMENT_EXPORT_COLUMNS: list[str] = [
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


EXPORT_COLUMN_RENAMES: Mapping[str, str] = {
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


EXPORT_COALESCE_SOURCES: Mapping[str, Sequence[str]] = {
    "OpenAlex.PMID": ["OpenAlex.PMID", "PubMed.PMID", "scholar.PMID"],
    "OpenAlex.DOI": ["OpenAlex.DOI", "PubMed.DOI", "scholar.DOI", "doi_normalised"],
    "crossref.DOI": ["crossref.DOI", "doi_normalised", "PubMed.DOI", "scholar.DOI"],
}


EXPORT_SORT_FALLBACK: list[str] = [
    "ChEMBL.document_chembl_id",
    "PubMed.PMID",
    "scholar.PMID",
    "OpenAlex.PMID",
    "ChEMBL.pubmed_id",
]


EXPORT_STREAM_CHUNK_SIZE = 10_000


def coalesce_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.Series:
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


def prepare_document_export_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Rename and project columns to match the export schema."""

    frame = df.copy()
    rename_map = {k: v for k, v in EXPORT_COLUMN_RENAMES.items() if k in frame.columns}
    if rename_map:
        frame = frame.rename(columns=rename_map)

    for target, sources in EXPORT_COALESCE_SOURCES.items():
        frame[target] = coalesce_columns(frame, sources)

    for column in DOCUMENT_EXPORT_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""

    return frame[DOCUMENT_EXPORT_COLUMNS]


def iter_document_export_chunks(
    df: pd.DataFrame, *, chunk_size: int | None
) -> Iterable[pd.DataFrame]:
    """Yield export-ready chunks with deterministic column ordering."""

    total_rows = len(df)
    if total_rows == 0:
        empty = dataframe_to_strings(df.copy(), skip=NUMERIC_EXPORT_COLUMNS)
        yield prepare_document_export_frame(empty)
        return

    effective_size = chunk_size if chunk_size and chunk_size > 0 else total_rows
    for start in range(0, total_rows, effective_size):
        stop = start + effective_size
        chunk = df.iloc[start:stop].copy()
        chunk = dataframe_to_strings(chunk, skip=NUMERIC_EXPORT_COLUMNS)
        yield prepare_document_export_frame(chunk)


def resolve_export_chunk_size(value: int | None) -> int | None:
    """Return ``value`` when positive, otherwise ``None``."""

    if value is None:
        return None
    if value <= 0:
        logger.warning("invalid_csv_chunksize", value=value)
        return None
    return value


def write_document_export_chunks(
    chunks: Iterable[pd.DataFrame],
    path: Path,
    *,
    cfg: Config,
    key_cols: Sequence[str],
    chunk_size: int | None,
) -> Path:
    """Stream ``chunks`` to ``path`` using the deterministic CSV writer."""

    if chunk_size:
        return write_csv_chunks_deterministic(
            chunks,
            path,
            col_order=list(DOCUMENT_EXPORT_COLUMNS),
            key_cols=list(key_cols),
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
            chunksize=chunk_size,
            merge_chunksize=chunk_size,
        )

    return write_csv_chunks_deterministic(
        chunks,
        path,
        col_order=list(DOCUMENT_EXPORT_COLUMNS),
        key_cols=list(key_cols),
        sep=cfg.io.csv_sep,
        encoding=cfg.io.csv_encoding,
        cfg=cfg,
    )


def load_export_ready_frame(path: Path, *, cfg: Config) -> pd.DataFrame:
    """Load the streamed CSV export for quality reporting."""

    try:
        loaded = pd.read_csv(
            path,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
    except pd.errors.EmptyDataError:
        loaded = pd.DataFrame(columns=list(DOCUMENT_EXPORT_COLUMNS))
    missing = [column for column in DOCUMENT_EXPORT_COLUMNS if column not in loaded.columns]
    for column in missing:
        loaded[column] = ""
    return loaded[DOCUMENT_EXPORT_COLUMNS]


def finalise_document_export(
    df: pd.DataFrame | Iterable[pd.DataFrame],
    output: Path,
    cfg: Config,
    *,
    input_csv: Path,
    key_columns: Sequence[str] | None = None,
    chunk_size: int | None = None,
) -> int:
    """Validate input frames and write CSV/metadata artefacts."""

    from schemas import DocumentsSchema

    if isinstance(df, pd.DataFrame):
        frames_iterable: Iterable[pd.DataFrame] = (df,)
    else:
        frames_iterable = df

    frames_iterator = iter(frames_iterable)
    analysis_iter, process_iter = tee(frames_iterator)

    required_cols = {
        name for name, col in DocumentsSchema.columns.items() if col.required
    }
    optional_cols = set(DocumentsSchema.columns) - required_cols

    present_columns: set[str] = set()
    for frame in analysis_iter:
        prepared = build_dataframe(
            add_pipeline_metadata(frame),
            columns=DOCUMENT_SCHEMA_COLUMNS,
            fill_missing=False,
        )
        present_columns.update(prepared.columns)

    missing_required = required_cols - present_columns
    missing_optional = optional_cols - present_columns

    if missing_required:
        logger.warning(
            "validation_skipped_missing_required",
            columns=sorted(missing_required),
        )
    elif missing_optional:
        logger.warning(
            "missing_optional_columns",
            columns=sorted(missing_optional),
        )

    should_validate = not missing_required

    stream_chunk = max(1, int(chunk_size or EXPORT_STREAM_CHUNK_SIZE))
    failure_path = output.with_name(f"{output.stem}_failure_cases.csv")
    errors = SidecarErrors()
    rows_total = 0
    rows_kept = 0
    exit_code = 0
    emitted_chunk = False

    def _validated_chunks() -> Iterator[pd.DataFrame]:
        nonlocal rows_total, rows_kept, exit_code, emitted_chunk
        for frame in process_iter:
            with_metadata = add_pipeline_metadata(frame)
            ordered = build_dataframe(
                with_metadata, columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False
            )
            rows_total += len(ordered)
            validated = ordered
            if should_validate:
                try:
                    validated = DocumentsSchema.validate(ordered, lazy=True)
                except SchemaErrors as exc:
                    for row in exc.failure_cases.to_dict("records"):
                        errors.add_error(row)
                    logger.error(
                        "document_validation_failed",
                        failure_count=len(exc.failure_cases),
                        failure_path=str(failure_path),
                        error=str(exc),
                    )
                    validated = getattr(exc, "validated_data", ordered)
                    exit_code = 1
            rows_kept += len(validated)
            cleaned = build_dataframe(
                validated, columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False
            )
            for chunk in iter_document_export_chunks(cleaned, chunk_size=stream_chunk):
                emitted_chunk = True
                yield chunk

        if not emitted_chunk:
            empty = build_dataframe(
                add_pipeline_metadata(pd.DataFrame()),
                columns=DOCUMENT_SCHEMA_COLUMNS,
                fill_missing=False,
            )
            for chunk in iter_document_export_chunks(empty, chunk_size=stream_chunk):
                emitted_chunk = True
                yield chunk

    export_generator = _validated_chunks()

    key_cols: list[str] = []
    if key_columns:
        for column in key_columns:
            mapped = EXPORT_COLUMN_RENAMES.get(column, column)
            if mapped in DOCUMENT_EXPORT_COLUMNS and mapped not in key_cols:
                key_cols.append(mapped)
    if not key_cols:
        for candidate in EXPORT_SORT_FALLBACK:
            if candidate in DOCUMENT_EXPORT_COLUMNS:
                key_cols = [candidate]
                break
    if not key_cols:
        key_cols = [DOCUMENT_EXPORT_COLUMNS[0]]

    col_order = list(DOCUMENT_EXPORT_COLUMNS)
    try:
        csv_path = write_csv_chunks_deterministic(
            export_generator,
            output,
            cfg=cfg,
            key_cols=key_cols,
            col_order=col_order,
            chunksize=stream_chunk,
            merge_chunksize=stream_chunk,
            sort_chunksize=stream_chunk,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
    except OSError as exc:
        logger.error("csv_write_failed", error=str(exc), path=str(output))
        return 1

    errors.save(failure_path)

    rows_dropped = rows_total - rows_kept
    logger.info("write_done", rows=rows_kept, path=str(csv_path))

    try:
        export_ready = load_export_ready_frame(csv_path, cfg=cfg)
    except (OSError, ValueError) as exc:
        logger.error("quality_frame_load_failed", error=str(exc), path=str(csv_path))
        return 1

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
        report = build_quality_report(export_ready)
        save_quality_report(report, quality_path)
    except (OSError, TypeError, ValueError) as exc:
        logger.error(
            "quality_report_write_failed",
            error=str(exc),
            path=str(quality_path),
        )
        return 1

    try:
        analyze_table_quality(export_ready, table_name=str(csv_path.with_suffix("")))
    except ValueError as exc:
        logger.error("quality_report_generation_failed", error=str(exc))
        return 1
    return exit_code


# ---------------------------------------------------------------------------
# PubMed aggregation


def build_fallback_doi_map(
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
        raise ValueError(f"missing columns in fallback DOI file: {', '.join(missing)}")

    mapping: dict[str, str] = {}
    for pmid_value, doi_value in frame[[pmid_column, doi_column]].itertuples(
        index=False
    ):
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
    return_generator: bool = False,
) -> pd.DataFrame | Iterator[pd.DataFrame]:
    """Retrieve metadata for a sequence of PubMed identifiers."""

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

    documents_limiter = get_limiter(
        "documents_global", rate_cfg.global_rps, rate_cfg.global_burst
    )

    def _service_limiter(
        name: str,
        *,
        rps: int | None,
        burst: int | None,
    ) -> RateLimiter | None:
        effective_rps = rps if rps is not None else rate_cfg.global_rps
        effective_burst = burst if burst is not None else rate_cfg.global_burst
        if (
            rps is None
            and burst is None
            or (
                effective_rps == rate_cfg.global_rps
                and effective_burst == rate_cfg.global_burst
            )
        ):
            return None
        return get_limiter(f"documents_{name}", effective_rps, effective_burst)

    pubmed_service_limiter = _service_limiter(
        "pubmed", rps=getattr(pubmed_cfg, "rps", None), burst=getattr(pubmed_cfg, "burst", None)
    )
    semantic_service_limiter = _service_limiter(
        "semantic_scholar",
        rps=getattr(semantic_scholar_cfg, "rps", None),
        burst=getattr(semantic_scholar_cfg, "burst", None),
    )
    openalex_service_limiter = _service_limiter(
        "openalex",
        rps=getattr(openalex_cfg, "rps", None),
        burst=getattr(openalex_cfg, "burst", None),
    )
    crossref_service_limiter = _service_limiter(
        "crossref",
        rps=getattr(crossref_cfg, "rps", None),
        burst=getattr(crossref_cfg, "burst", None),
    )

    def _acquire_documents(limiter: RateLimiter | None) -> None:
        documents_limiter.acquire()
        if limiter is not None:
            limiter.acquire()

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

    def _executor_capacity(limiter: RateLimiter | None, burst: int | None) -> int:
        limit = burst if burst is not None else 1
        if limiter is not None:
            limit = min(limit, limiter.burst)
        return max(1, limit)

    def _fetch_batch(
        first: object,
        *rest: object,
        __cfg: Config = session_cfg,
        **__: object,
    ) -> list[dict[str, str]]:
        """Fetch metadata for a batch of PMIDs."""

        batch_list = _coerce_batch_argument(first, *rest)

        def _summarise_batch(pmids: Sequence[str]) -> dict[str, object]:
            sample_limit = 5
            sample = [pmids[index] for index in range(min(len(pmids), sample_limit))]
            if len(pmids) > sample_limit:
                sample.append("...")
            return {"pmids_count": len(pmids), "pmids_sample": sample}

        batch_summary = _summarise_batch(batch_list)

        def _invoke_with_session(
            fetcher: Callable[
                [requests.Session, str, Any, RateLimiter | None], dict[str, str]
            ],
            identifier: str,
            *,
            cfg_obj: Any,
            limiter: RateLimiter | None,
        ) -> dict[str, str]:
            with session_with_retry(__cfg.api, __cfg.retry) as nested_session:
                return fetcher(nested_session, identifier, cfg_obj, limiter)

        try:
            with session_with_retry(__cfg.api, __cfg.retry) as session:
                _acquire_documents(pubmed_service_limiter)
                pubmed_list = pl.fetch_pubmed_batch(
                    session, batch_list, sleep, cfg=pubmed_cfg
                )

                pmids_in_batch = [p.get("PubMed.PMID", "") for p in pubmed_list]
                semantic_pmids = [pmid for pmid in pmids_in_batch if pmid]

                semsch_map: dict[str, dict[str, str]] = {}
                if semantic_pmids:
                    _acquire_documents(semantic_service_limiter)
                    semsch_list = ssl.fetch_semantic_scholar_batch(
                        session, semantic_pmids, sleep, cfg=semantic_scholar_cfg
                    )

                    semsch_map = {
                        s.get("scholar.PMID"): s
                        for s in semsch_list
                        if s.get("scholar.PMID")
                    }

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
                        _acquire_documents(semantic_service_limiter)
                        fallback_record = ssl.fetch_semantic_scholar(
                            session, pmid, sleep, cfg=semantic_scholar_cfg
                        )
                        semsch_map[pmid] = fallback_record

                combined_records: list[dict[str, str]] = []

                plan: list[tuple[int, dict[str, str], dict[str, str], str, str]] = []
                openalex_lookup: dict[str, list[int]] = {}
                crossref_lookup: dict[str, list[int]] = {}
                openalex_total = 0
                crossref_total = 0

                for index, pubmed in enumerate(pubmed_list):
                    pmid = pubmed.get("PubMed.PMID", "")
                    semsch = semsch_map.get(pmid, {}) if pmid else {}
                    fallback_doi = ""
                    if fallback_doi_map:
                        fallback_doi = fallback_doi_map.get(pmid, "")
                    doi = (
                        pubmed.get("PubMed.DOI")
                        or semsch.get("scholar.DOI")
                        or fallback_doi
                        or ""
                    )
                    plan.append((index, pubmed, semsch, pmid, doi))
                    if pmid:
                        openalex_lookup.setdefault(pmid, []).append(index)
                        openalex_total += 1
                    crossref_key = doi or ""
                    crossref_lookup.setdefault(crossref_key, []).append(index)
                    crossref_total += 1

                openalex_results: dict[int, dict[str, str]] = {}
                crossref_results: dict[int, dict[str, str]] = {}

                openalex_jobs = list(openalex_lookup.keys())

                def _fetch_openalex_job(pmid: str) -> dict[str, str]:
                    _acquire_documents(openalex_service_limiter)
                    return _invoke_with_session(
                        ocl.fetch_openalex,
                        pmid,
                        cfg_obj=openalex_cfg,
                        limiter=openalex_limiter,
                    )

                if openalex_total:
                    logger.info(
                        "documents_cache_reuse",
                        service="openalex",
                        total=openalex_total,
                        unique=len(openalex_jobs),
                        hits=openalex_total - len(openalex_jobs),
                    )

                openalex_by_key: dict[str, dict[str, str]] = {}

                if openalex_jobs:
                    if openalex_executor is None:
                        for pmid in openalex_jobs:
                            openalex_by_key[pmid] = _fetch_openalex_job(pmid)
                    else:
                        future_to_key = {
                            openalex_executor.submit(_fetch_openalex_job, pmid): pmid
                            for pmid in openalex_jobs
                        }
                        for future in as_completed(future_to_key):
                            key = future_to_key[future]
                            openalex_by_key[key] = future.result()

                for pmid, indexes in openalex_lookup.items():
                    result = openalex_by_key.get(pmid, {})
                    for idx in indexes:
                        openalex_results[idx] = result

                crossref_jobs = list(crossref_lookup.keys())

                def _fetch_crossref_job(doi: str) -> dict[str, str]:
                    _acquire_documents(crossref_service_limiter)
                    return _invoke_with_session(
                        ocl.fetch_crossref,
                        doi,
                        cfg_obj=crossref_cfg,
                        limiter=crossref_limiter,
                    )

                if crossref_total:
                    logger.info(
                        "documents_cache_reuse",
                        service="crossref",
                        total=crossref_total,
                        unique=len(crossref_jobs),
                        hits=crossref_total - len(crossref_jobs),
                    )

                crossref_by_key: dict[str, dict[str, str]] = {}

                if crossref_jobs:
                    if crossref_executor is None:
                        for doi in crossref_jobs:
                            crossref_by_key[doi] = _fetch_crossref_job(doi)
                    else:
                        future_to_key = {
                            crossref_executor.submit(_fetch_crossref_job, doi): doi
                            for doi in crossref_jobs
                        }
                        for future in as_completed(future_to_key):
                            key = future_to_key[future]
                            crossref_by_key[key] = future.result()

                for doi, indexes in crossref_lookup.items():
                    result = crossref_by_key.get(doi, {})
                    for idx in indexes:
                        crossref_results[idx] = result

                for index, pubmed, semsch, pmid, _ in plan:
                    openalex = openalex_results.get(index, {}) if pmid else {}
                    crossref = crossref_results.get(index, {})
                    combined = merge_metadata(pubmed, semsch, openalex, crossref)
                    combined_records.append(combined)
                return combined_records
        except requests.RequestException as exc:  # pragma: no cover - network errors
            logger.warning(
                "pubmed_batch_request_failed",
                **batch_summary,
                error=str(exc),
            )
            return _failure_records(batch_list, str(exc))
        except Exception as exc:  # pragma: no cover - defensive safety net
            logger.warning(
                "pubmed_batch_unexpected_error",
                **batch_summary,
                error=str(exc),
            )
            return _failure_records(batch_list, str(exc))

    openalex_capacity = _executor_capacity(openalex_limiter, openalex_cfg.burst)
    crossref_capacity = _executor_capacity(crossref_limiter, crossref_cfg.burst)

    openalex_executor: ThreadPoolExecutor | None = None
    crossref_executor: ThreadPoolExecutor | None = None

    iterator = (p for p in pmids if p)

    tasks: dict[Future[list[dict[str, str]]], tuple[int, list[str]]] = {}
    completed: dict[int, list[dict[str, str]]] = {}
    next_to_emit = 0
    processed = 0
    max_in_flight = max(1, max_workers * 2)

    with ExitStack() as stack:
        batch_executor = stack.enter_context(
            ThreadPoolExecutor(max_workers=max_workers)
        )
        if openalex_capacity > 1:
            openalex_executor = stack.enter_context(
                ThreadPoolExecutor(max_workers=openalex_capacity)
            )
        if crossref_capacity > 1:
            crossref_executor = stack.enter_context(
                ThreadPoolExecutor(max_workers=crossref_capacity)
            )

        offset = 0
        pending: set[Future[list[dict[str, str]]]] = set()

        def _drain_future(
            done_future: Future[list[dict[str, str]]],
        ) -> Iterator[list[dict[str, str]]]:
            nonlocal processed, next_to_emit

            pending.remove(done_future)
            batch_id, batch_pmids = tasks.pop(done_future)
            completed[batch_id] = done_future.result()
            processed += len(batch_pmids)
            logger.info("documents_processed", count=processed)

            while next_to_emit in completed:
                yield completed.pop(next_to_emit)
                next_to_emit += 1

        def _emit_ready_batches() -> Iterator[list[dict[str, str]]]:
            nonlocal next_to_emit

            while next_to_emit in completed:
                yield completed.pop(next_to_emit)
                next_to_emit += 1

        def _iter_records() -> Iterator[list[dict[str, str]]]:
            nonlocal offset

            for batch in _chunked(iterator, batch_size):
                if not batch:
                    continue
                future = batch_executor.submit(_fetch_batch, batch)
                tasks[future] = (offset, batch)
                pending.add(future)
                offset += len(batch)
                if len(pending) >= max_in_flight:

                    done_future = next(as_completed(list(pending)))
                    yield from _drain_future(done_future)

                yield from _emit_ready_batches()

            for done_future in as_completed(pending):
                yield from _drain_future(done_future)

            pending.clear()

            yield from _emit_ready_batches()

        def _iter_frames() -> Iterator[pd.DataFrame]:
            for records_batch in _iter_records():
                if not records_batch:
                    yield build_dataframe([], columns=DOCUMENT_SCHEMA_COLUMNS)
                    continue
                yield build_dataframe(
                    records_batch, columns=DOCUMENT_SCHEMA_COLUMNS
                )

    frame_iter = _iter_frames()
    if return_generator:
        return frame_iter

    frame_iter, concat_iter = tee(frame_iter)
    try:
        first_frame = next(concat_iter)
    except StopIteration:
        return build_dataframe([], columns=DOCUMENT_SCHEMA_COLUMNS)
    return pd.concat(chain([first_frame], concat_iter), ignore_index=True)


def run_pubmed_pipeline(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``pubmed`` sub-command using reusable helpers."""

    from schemas import normalize_documents

    pubmed_defaults = cfg.document.pubmed
    limit = getattr(args, "limit", pubmed_defaults.limit)
    if limit is not None and limit < 0:
        logger.error(
            "invalid_limit",
            section="document.pubmed",
            limit=limit,
        )
        return 1
    try:
        pmids_iter = io.read_ids(
            args.input_csv,
            column=getattr(args, "column", pubmed_defaults.column),
            cfg=cfg.io,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "input_read_failed",
            error=str(exc),
            path=str(args.input_csv),
        )
        return 1
    offset = getattr(args, "offset", 0)
    if offset:
        pmids_iter = islice(pmids_iter, offset, None)
        logger.info("process_offset", offset=offset)
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
            fallback_map = build_fallback_doi_map(
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
        frame_iter = fetch_pubmed_records(
            pmids,
            cfg,
            sleep=getattr(args, "sleep", pubmed_defaults.sleep),
            pubmed_cfg=cfg.pubmed,
            semantic_scholar_cfg=cfg.semantic_scholar,
            openalex_cfg=cfg.openalex,
            crossref_cfg=cfg.crossref,
            max_workers=getattr(args, "workers", pubmed_defaults.workers),
            batch_size=getattr(args, "batch_size", pubmed_defaults.batch_size),
            fallback_doi_map=fallback_doi_map,
            return_generator=True,
        )
        output = Path(args.output_csv or io.default_output_path(args.input_csv, cfg.io))
        normalised_frames = (normalize_documents(frame) for frame in frame_iter)
        exit_code = finalise_document_export(
            normalised_frames,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
            chunk_size=getattr(args, "batch_size", pubmed_defaults.batch_size),
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("pubmed_pipeline_failed", error=str(exc))
        return 1
    return exit_code


def run_chembl_pipeline(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``chembl`` sub-command using reusable helpers."""

    from schemas import normalize_documents

    chembl_defaults = cfg.document.chembl
    limit = getattr(args, "limit", chembl_defaults.limit)
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="document.chembl", limit=limit)
        return 1

    with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
        try:
            ids_iter = io.read_ids(
                args.input_csv,
                column=getattr(args, "column", chembl_defaults.column),
                cfg=cfg.io,
            )
        except (FileNotFoundError, ValueError) as exc:
            logger.error(
                "input_read_failed",
                error=str(exc),
                path=str(args.input_csv),
            )
            return 1

        offset = getattr(args, "offset", 0)
        if offset:
            ids_iter = islice(ids_iter, offset, None)
            logger.info("process_offset", offset=offset)

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
                chunk_size=getattr(args, "chunk_size", chembl_defaults.chunk_size),
                timeout=getattr(args, "timeout", chembl_defaults.timeout),
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error(
                "chembl_documents_fetch_failed",
                error=str(exc),
                chunk_size=getattr(args, "chunk_size", chembl_defaults.chunk_size),
            )
            return 1
        if "doi" in df.columns:
            df["doi"] = df["doi"].map(normalise_doi)
        output = Path(args.output_csv or io.default_output_path(args.input_csv, cfg.io))
        df = normalize_documents(df)
        exit_code = finalise_document_export(
            df,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
            chunk_size=getattr(args, "chunk_size", chembl_defaults.chunk_size),
        )
        return exit_code


def run_all_pipeline(cfg: Config, args: argparse.Namespace) -> int:
    """Run ChEMBL and PubMed pipelines and merge their outputs."""

    from schemas import normalize_documents

    all_defaults = cfg.document.all
    limit = getattr(args, "limit", all_defaults.limit)
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="document.all", limit=limit)
        return 1

    try:
        ids_iter = io.read_ids(
            args.input_csv,
            column=getattr(args, "column", all_defaults.column),
            cfg=cfg.io,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "input_read_failed",
            error=str(exc),
            path=str(args.input_csv),
        )
        return 1

    offset = getattr(args, "offset", 0)
    if offset:
        ids_iter = islice(ids_iter, offset, None)
        logger.info("process_offset", offset=offset)

    ids_source: Iterable[str] = ids_iter
    if limit is not None:
        limited_ids = list(islice(ids_iter, limit))
        ids_source = limited_ids
        logger.info("process_limit", limit=len(limited_ids))

    iterator = iter(ids_source)
    sample_size = getattr(args, "chunk_size", all_defaults.chunk_size)
    sample_ids = list(islice(iterator, sample_size))
    ids_for_fetch = chain(sample_ids, iterator)
    try:
        with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
            doc_df = cl.get_documents(
                ids_for_fetch,
                cfg=cfg.api,
                client=client,
                chunk_size=getattr(args, "chunk_size", all_defaults.chunk_size),
                timeout=getattr(args, "timeout", all_defaults.timeout),
            )
    except (requests.RequestException, ValueError) as exc:
        logger.error(
            "chembl_documents_fetch_failed",
            ids=sample_ids,
            error=str(exc),
            chunk_size=getattr(args, "chunk_size", all_defaults.chunk_size),
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
        exit_code = finalise_document_export(
            processed,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
            chunk_size=getattr(args, "chunk_size", all_defaults.chunk_size),
        )
        return exit_code

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
    pubmed_frames = fetch_pubmed_records(
        pmids,
        cfg,
        sleep=getattr(args, "sleep", all_defaults.sleep),
        semantic_scholar_cfg=cfg.semantic_scholar,
        openalex_cfg=cfg.openalex,
        crossref_cfg=cfg.crossref,
        max_workers=getattr(args, "workers", all_defaults.workers),
        batch_size=getattr(args, "batch_size", all_defaults.batch_size),
        pubmed_cfg=cfg.pubmed,
        fallback_doi_map=doi_fallback_map or None,
        return_generator=True,
    )
    concat_iter = iter(pubmed_frames)
    try:
        first_frame = next(concat_iter)
    except StopIteration:
        pub_df = build_dataframe([], columns=DOCUMENT_SCHEMA_COLUMNS)
    else:
        pub_df = pd.concat(chain([first_frame], concat_iter), ignore_index=True)
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
    exit_code = finalise_document_export(
        processed,
        output,
        cfg,
        input_csv=Path(args.input_csv),
        key_columns=["document_chembl_id"],
        chunk_size=getattr(args, "chunk_size", all_defaults.chunk_size),
    )
    return exit_code


__all__ = [
    "CH_EMBL_COLUMNS",
    "PUBMED_COLUMNS",
    "SEMANTIC_SCHOLAR_COLUMNS",
    "OPENALEX_COLUMNS",
    "CROSSREF_COLUMNS",
    "DOCUMENT_SCHEMA_COLUMNS",
    "merge_metadata",
    "build_dataframe",
    "merge_with_chembl",
    "dataframe_to_strings",
    "normalise_doi",
    "build_quality_report",
    "save_quality_report",
    "NUMERIC_EXPORT_COLUMNS",
    "DOCUMENT_EXPORT_COLUMNS",
    "EXPORT_COLUMN_RENAMES",
    "EXPORT_COALESCE_SOURCES",
    "EXPORT_SORT_FALLBACK",
    "EXPORT_STREAM_CHUNK_SIZE",
    "coalesce_columns",
    "prepare_document_export_frame",
    "iter_document_export_chunks",
    "resolve_export_chunk_size",
    "write_document_export_chunks",
    "load_export_ready_frame",
    "finalise_document_export",
    "build_fallback_doi_map",
    "fetch_pubmed_records",
    "run_pubmed_pipeline",
    "run_chembl_pipeline",
    "run_all_pipeline",
]
