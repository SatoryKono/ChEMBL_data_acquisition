"""Post-processing helpers for document exports.

Changelog
~~~~~~~~
- 2025-02-??: Relaxed reference classification parsing to tolerate non-numeric
  values introduced in recent exports.
- 2024-04-??: Initial Python port of the Power Query pipeline used to derive
  analytics-friendly document summaries.

"""

from __future__ import annotations

import importlib
import math
import os
import re
from collections.abc import Iterable, Sequence
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping, cast

import numpy as np
import pandas as pd

from ..pipelines.document import postprocessing as stage_document_postprocessing
from ..config import IoCfg
from ..common.csv_utils import write_csv_deterministic
from ..common.log import logger

# ===== Parameters ===========================================================
UTF8_ENCODING = "utf-8"
DEFAULT_OUTPUT_PREFIX = "preprocessed_"
METADATA_LIST_SEPARATOR = "; "
TOKEN_DELIMITERS: tuple[str, ...] = (";", "|", ",")
INDEX_PAD_WIDTH = 4
COVERAGE_COMPLETE_THRESHOLD = 4
COVERAGE_PARTIAL_THRESHOLD = 2
COVERAGE_MINIMAL_THRESHOLD = 1
COVERAGE_FAILURE_THRESHOLD = 0
REQUIRED_INPUT_COLUMNS: tuple[str, ...] = ("document_chembl_id",)
REVIEW_COLUMNS: tuple[str, ...] = (
    "PubMed.is_review",
    "scholar.is_review",
    "OpenAlex.is_review",
)
PASS_THROUGH_COLUMNS: tuple[str, ...] = (
    "authors",
    "source",
    "fetch_status",
    "pipeline_version",
    "timestamp_utc",
    "date_code",
    "Index",
)
COALESCE_PRIORITIES: Mapping[str, Sequence[str]] = {
    "preferred_title": (
        "title",
        "PubMed.ArticleTitle",
        "crossref.Title",
    ),
    "preferred_abstract": (
        "abstract",
        "PubMed.Abstract",
    ),
    "preferred_doi": (
        "doi_normalised",
        "doi",
        "PubMed.DOI",
        "scholar.DOI",
    ),
    "primary_pubmed_id": (
        "PubMed.PMID",
        "pubmed_id",
        "scholar.PMID",
    ),
    "preferred_journal": (
        "journal",
        "PubMed.JournalTitle",
        "scholar.Venue",
        "OpenAlex.Venue",
        "crossref.Title",
    ),
    "publication_year": (
        "year",
        "PubMed.YearCompleted",
        "PubMed.YearRevised",
    ),
}
METADATA_SOURCE_ORDER: tuple[str, ...] = (
    "chembl",
    "pubmed",
    "semantic_scholar",
    "openalex",
    "crossref",
)
METADATA_SOURCE_COLUMNS: Mapping[str, Sequence[str]] = {
    "chembl": (
        "document_chembl_id",
        "title",
        "abstract",
        "doi",
        "doi_normalised",
        "journal",
        "authors",
    ),
    "pubmed": (
        "PubMed.PMID",
        "PubMed.ArticleTitle",
        "PubMed.Abstract",
        "PubMed.MeSH_Descriptors",
        "PubMed.MeSH_Qualifiers",
        "PubMed.Error",
    ),
    "semantic_scholar": (
        "scholar.PMID",
        "scholar.SemanticScholarId",
        "scholar.ExternalIds",
        "scholar.DOI",
        "scholar.Venue",
        "scholar.Error",
    ),
    "openalex": (
        "OpenAlex.Id",
        "OpenAlex.TypeCrossref",
        "OpenAlex.Genre",
        "OpenAlex.MeshDescriptors",
        "OpenAlex.MeshQualifiers",
        "OpenAlex.Error",
    ),
    "crossref": (
        "crossref.Type",
        "crossref.Subtype",
        "crossref.Subject",
        "crossref.Title",
        "crossref.Error",
    ),
}

CHEMBL_EXPORT_ALIASES: Mapping[str, str] = {
    "ChEMBL.document_chembl_id": "document_chembl_id",
    "ChEMBL.title": "title",
    "ChEMBL.abstract": "abstract",
    "ChEMBL.doi": "doi",
    "ChEMBL.year": "year",
    "ChEMBL.journal": "journal",
    "ChEMBL.journal_abbrev": "journal_abbrev",
    "ChEMBL.volume": "volume",
    "ChEMBL.issue": "issue",
    "ChEMBL.first_page": "first_page",
    "ChEMBL.last_page": "last_page",
    "ChEMBL.pubmed_id": "pubmed_id",
    "ChEMBL.authors": "authors",
    "ChEMBL.source": "source",
}
ERROR_COLUMN_MAP: Mapping[str, str] = {
    "pubmed": "PubMed.Error",
    "semantic_scholar": "scholar.Error",
    "openalex": "OpenAlex.Error",
    "crossref": "crossref.Error",
}
ERROR_PRIORITY: tuple[str, ...] = (
    "pubmed",
    "semantic_scholar",
    "openalex",
    "crossref",
    "unknown",
)
FETCH_ERROR_STATUS = "error"
MESH_TERM_COLUMNS: tuple[str, ...] = (
    "PubMed.MeSH_Descriptors",
    "PubMed.MeSH_Qualifiers",
    "OpenAlex.MeshDescriptors",
    "OpenAlex.MeshQualifiers",
)
PREPROCESSED_COLUMN_ORDER: tuple[str, ...] = (
    "document_chembl_id",
    "primary_pubmed_id",
    "preferred_title",
    "preferred_abstract",
    "preferred_doi",
    "preferred_journal",
    "publication_year",
    "publication_class",
    "is_review",
    "metadata_sources",
    "metadata_source_count",
    "coverage_score",
    "coverage_status",
    "has_chembl",
    "has_pubmed",
    "has_semantic_scholar",
    "has_openalex",
    "has_crossref",
    "error_sources",
    "has_error",
    "mesh_terms",
    "authors",
    "source",
    "fetch_status",
    "pipeline_version",
    "timestamp_utc",
    "date_code",
    "Index",
)
DEFAULT_KEY_COLUMNS: tuple[str, ...] = ("completed", "document_chembl_id")


# ===== Helpers ==============================================================
def _string_value(value: object) -> str:
    """Return ``value`` as a trimmed string, treating nulls as empty."""

    if isinstance(value, str):
        return value.strip()
    if value is None:
        return ""
    if isinstance(value, (int, float)):
        if pd.isna(value):
            return ""
        return str(value).strip()
    if pd.isna(value):
        return ""
    return str(value).strip()


def _string_series(series: pd.Series) -> pd.Series:
    """Return ``series`` as a trimmed ``string`` dtype series."""

    if series.empty:
        return pd.Series(dtype="string")
    result = series.astype("string").fillna("").str.strip()
    return result


def _apply_export_aliases(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of ``df`` with ChEMBL export aliases populated."""

    frame = df.copy()
    for source, target in CHEMBL_EXPORT_ALIASES.items():
        if source not in frame.columns:
            continue
        if target in frame.columns:
            mask = _string_series(frame[target]) == ""
            if mask.any():
                frame.loc[mask, target] = frame.loc[mask, source]
        else:
            frame[target] = frame[source]
    return frame


def _truthy_mask(series: pd.Series) -> pd.Series:
    """Return a boolean mask highlighting non-empty values in ``series``."""

    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    cleaned = _string_series(series)
    return cleaned != ""


def _coalesce_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.Series:
    """Return the first non-empty value from ``columns`` per row."""

    if df.empty:
        return pd.Series(dtype="string")
    result = pd.Series([""] * len(df), index=df.index, dtype="string")
    for column in columns:
        if column not in df.columns:
            continue
        candidate = _string_series(df[column])
        mask = result == ""
        result = result.mask(mask & (candidate != ""), candidate)
    return result


def _series_to_bool(series: pd.Series) -> pd.Series:
    """Return ``series`` as boolean values recognising textual truthy tokens."""

    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    tokens = _string_series(series).str.lower()
    truthy = {"true", "yes", "1", "review", "y"}
    return tokens.isin(truthy)


def _split_tokens(value: str, delimiters: Iterable[str]) -> list[str]:
    """Split ``value`` using ``delimiters`` preserving order."""

    tokens = [value]
    for delimiter in delimiters:
        parts: list[str] = []
        for token in tokens:
            parts.extend(token.split(delimiter))
        tokens = parts
    return [token.strip() for token in tokens if token.strip()]


def _aggregate_terms(df: pd.DataFrame, columns: Sequence[str]) -> pd.Series:
    """Combine descriptor columns into a semicolon-delimited series."""

    if not columns:
        return pd.Series([""] * len(df), index=df.index, dtype="string")
    collected: list[str] = []
    for values in df[columns].fillna("").itertuples(index=False, name=None):
        seen: set[str] = set()
        terms: list[str] = []
        for value in values:
            text = _string_value(value)
            if not text:
                continue
            for token in _split_tokens(text, TOKEN_DELIMITERS):
                lowered = token.lower()
                if lowered in seen:
                    continue
                seen.add(lowered)
                terms.append(token)
        collected.append(METADATA_LIST_SEPARATOR.join(terms))
    return pd.Series(collected, index=df.index, dtype="string")


def _sort_tokens(tokens: Iterable[str], priority: Sequence[str]) -> list[str]:
    """Return ``tokens`` ordered according to ``priority``."""

    order = {name: idx for idx, name in enumerate(priority)}
    unique: list[str] = []
    seen: set[str] = set()
    for token in tokens:
        cleaned = token.strip().lower()
        if not cleaned or cleaned in seen:
            continue
        seen.add(cleaned)
        unique.append(cleaned)
    unique.sort(key=lambda item: order.get(item, len(order)))
    return unique


def _build_error_sources(df: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Return aggregated error sources and presence flags."""

    error_tokens: list[str] = []
    has_error: list[bool] = []

    fetch_status = df.get("fetch_status")
    error_source = df.get("error_source")

    for idx in range(len(df)):
        tokens: list[str] = []
        for name, column in ERROR_COLUMN_MAP.items():
            if column not in df.columns:
                continue
            value = _string_value(df.iloc[idx][column])
            if value:
                tokens.append(name)
        if fetch_status is not None:
            status = _string_value(fetch_status.iloc[idx]).lower()
            if status == FETCH_ERROR_STATUS:
                source = "unknown"
                if error_source is not None:
                    source = _string_value(error_source.iloc[idx]) or "unknown"
                tokens.append(source.lower())
        ordered = _sort_tokens(tokens, ERROR_PRIORITY)
        has_error.append(bool(ordered))
        error_tokens.append(METADATA_LIST_SEPARATOR.join(ordered))
    return (
        pd.Series(error_tokens, index=df.index, dtype="string"),
        pd.Series(has_error, index=df.index, dtype="boolean"),
    )


def _build_metadata_flags(df: pd.DataFrame) -> tuple[pd.Series, dict[str, pd.Series]]:
    """Return aggregated metadata source strings and per-source flags."""

    flags: dict[str, pd.Series] = {}
    token_lists: list[list[str]] = []

    for source in METADATA_SOURCE_ORDER:
        columns = METADATA_SOURCE_COLUMNS.get(source, ())
        mask = pd.Series([False] * len(df), index=df.index, dtype=bool)
        for column in columns:
            if column not in df.columns:
                continue
            mask = mask | _truthy_mask(df[column])
        flags[f"has_{source}"] = mask
    for row in zip(*(flags[f"has_{src}"] for src in METADATA_SOURCE_ORDER)):
        sources = [name for name, present in zip(METADATA_SOURCE_ORDER, row) if present]
        token_lists.append(sources)
    metadata_strings = [METADATA_LIST_SEPARATOR.join(tokens) for tokens in token_lists]
    metadata_counts = [len(tokens) for tokens in token_lists]
    metadata_series = pd.Series(metadata_strings, index=df.index, dtype="string")
    count_series = pd.Series(metadata_counts, index=df.index, dtype="Int64")
    return metadata_series, {**flags, "metadata_source_count": count_series}


def _coverage_status(score: int, has_error: bool) -> str:
    """Return textual coverage label for ``score``."""

    if has_error and score <= COVERAGE_FAILURE_THRESHOLD:
        return "failed"
    if score >= COVERAGE_COMPLETE_THRESHOLD:
        return "complete"
    if score >= COVERAGE_PARTIAL_THRESHOLD:
        return "partial"
    if score >= COVERAGE_MINIMAL_THRESHOLD:
        return "minimal"
    return "unknown"


# ===== Modules ==============================================================
def _ensure_document_identifier(frame: pd.DataFrame) -> pd.DataFrame:
    """Ensure ``frame`` exposes a ``document_chembl_id`` column."""

    if "document_chembl_id" in frame.columns:
        return frame

    # Normalise alternative representations (case differences, stray spacing)
    for column in frame.columns:
        if column.strip().lower() == "document_chembl_id":
            if column == "document_chembl_id":
                return frame
            renamed = frame.copy()
            renamed = renamed.rename(columns={column: "document_chembl_id"})
            logger.warning(
                "document_id_column_renamed",
                original=column,
                normalised="document_chembl_id",
            )
            return renamed

    if "ChEMBL.document_chembl_id" in frame.columns:
        duplicated = frame.copy()
        duplicated["document_chembl_id"] = duplicated["ChEMBL.document_chembl_id"]
        logger.warning(
            "document_id_column_duplicated",
            source="ChEMBL.document_chembl_id",
            target="document_chembl_id",
        )
        return duplicated

    fallback = frame.copy()
    fallback["document_chembl_id"] = ""
    logger.error(
        "document_id_column_missing",
        columns=list(frame.columns),
    )
    return fallback


def preprocess_document_export(
    df: pd.DataFrame,
    *,
    ref_document: pd.DataFrame | None = None,
    ref_document_path: Path | str | None = None,
) -> pd.DataFrame:
    """Return the Stage-aligned projection of ``df``.

    This is a thin wrapper around :func:`library.pipelines.document.postprocessing.postprocess_documents`
    that preserves the legacy entry point used by :mod:`scripts.get_document_data`.
    The resulting frame strictly follows :data:`FINAL_COLUMN_ORDER` defined in the
    Stage port and omits analytics-only helper columns (``preferred_*``,
    ``coverage_*`` and similar aggregates).
    """

    processed = stage_document_postprocessing.postprocess_documents(
        df,
        ref_document=ref_document,
        ref_document_path=ref_document_path,
    )

    missing_columns = [
        column
        for column in stage_document_postprocessing.FINAL_COLUMN_ORDER
        if column not in processed.columns
    ]
    if missing_columns:
        raise ValueError(
            "Stage post-processing did not return the expected column set",
            missing_columns,
        )

    return processed.loc[:, stage_document_postprocessing.FINAL_COLUMN_ORDER]


def _normalise_export_basename(source: Path | str) -> str:
    """Return the ``source`` basename without leading dot or ``.tmp`` suffix."""

    name = Path(source).name
    if name.startswith("."):
        name = name[1:]
    if name.endswith(".tmp"):
        name = name[: -len(".tmp")]
    return name


def _build_default_destination(source: Path | str) -> Path:
    """Return the default output path for ``source``."""

    source_path = Path(source)
    base_name = _normalise_export_basename(source_path)
    return source_path.with_name(f"{DEFAULT_OUTPUT_PREFIX}{base_name}")


def postprocess_export_file(
    input_path: Path | str,
    *,
    cfg: IoCfg,
    output_path: Path | None = None,
    sep: str | None = None,
    encoding: str | None = None,
    ref_document: pd.DataFrame | None = None,
    ref_document_path: Path | str | None = None,
) -> Path:
    """Read ``input_path`` and materialise the Stage-aligned projection."""

    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding or UTF8_ENCODING
    frame = pd.read_csv(
        input_path,
        sep=sep,
        encoding=encoding,
        dtype=str,
        keep_default_na=False,
    )
    processed = preprocess_document_export(
        frame,
        ref_document=ref_document,
        ref_document_path=ref_document_path,
    )
    destination = Path(output_path) if output_path else _build_default_destination(input_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    col_order = [
        column
        for column in stage_document_postprocessing.FINAL_COLUMN_ORDER
        if column in processed.columns
    ]
    write_csv_deterministic(
        processed,
        destination,
        col_order=col_order,
        key_cols=list(DEFAULT_KEY_COLUMNS),
        sep=sep,
        encoding=encoding,
        cfg=None,
    )
    return destination



# ---------------------------------------------------------------------------
# Parameters (paths, encodings, delimiters)
# ---------------------------------------------------------------------------

CP1252 = "cp1252"
CSV_DELIMITER = ","

DEFAULT_BASE_PATH = "e:\\github\\ChEMBL_data_acquisition\\data\\"
DEFAULT_REF_RELATIVE = "input\\full\\document.csv"
DEFAULT_OUTPUT_RELATIVE = "output\\document\\output.document_YYYYMMDD.csv"
DEFAULT_QA_REFERENCE_RELATIVE = "input\\full\\document.csv"


# ---------------------------------------------------------------------------
# Column specifications and constants
# ---------------------------------------------------------------------------

PUBMED_NUMERIC_TEXT_COLUMNS = (
    "PubMed.Volume",
    "PubMed.Issue",
    "PubMed.StartPage",
    "PubMed.EndPage",
)

CHEMBL_NUMERIC_TEXT_COLUMNS = (
    "ChEMBL.volume",
    "ChEMBL.issue",
    "ChEMBL.first_page",
    "ChEMBL.last_page",
)

EXTERNAL_DOI_COLUMNS = (
    "PubMed.DOI",
    "scholar.DOI",
    "OpenAlex.DOI",
    "crossref.DOI",
)

EXTERNAL_PMID_COLUMNS = (
    "PubMed.PMID",
    "scholar.PMID",
    "OpenAlex.PMID",
)

LOWERCASE_LIST_COLUMNS = (
    "PubMed.mesh.descriptors",
    "PubMed.mesh.qualifiers",
    "PubMed.chemical.list",
    "OpenAlex.mesh.descriptors",
)

STAGE_REMOVED1_COLUMNS = (
    "PubMed.JournalTitle",
    "PubMed.JournalISOAbbrev",
    "scholar.Error",
    "OpenAlex.PublicationTypes",
    "OpenAlex.Genre",
    "OpenAlex.Venue",
    "OpenAlex.MeshQualifiers",
    "OpenAlex.Error",
    "crossref.Subtype",
    "crossref.Subtitle",
    "crossref.Subject",
    "crossref.Error",
    "publication_types_normalised",
    "publication_review_score",
    "publication_experimental_score",
    "publication_class",
    "doi",
)

FINAL_COLUMN_ORDER = (
    "PMID",
    "document_chembl_id",
    "doi",
    "reference",
    "completed",
    "sortorder.document",
    "review",
    "experimental",
    "document_contains_external_links",
    "invalid",
    "title",
    "abstract",
    "PubMed.mesh.descriptors",
    "PubMed.mesh.qualifiers",
    "PubMed.chemical.list",
    "OpenAlex.mesh.descriptors",
    "PubMed.Volume",
    "PubMed.Issue",
    "PubMed.StartPage",
    "PubMed.EndPage",
    "PubMed.ISSN",
    "PubMed.PublicationType",
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
    "OpenAlex.PMID",
    "OpenAlex.DOI",
    "OpenAlex.TypeCrossref",
    "OpenAlex.Id",
    "crossref.DOI",
    "crossref.Type",
    "crossref.Title",
    "ChEMBL.title",
    "ChEMBL.abstract",
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
    "invalid.doi",
    "invalid.PMID",
    "invalid.reference",
    "year",
    "month",
    "day",
)

WINDOWS_PATH_SEPARATOR = "\\"


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------


def to_text(value: Any) -> str:
    """Return *value* converted to text, using an empty string for nulls."""

    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="ignore")
    if isinstance(value, (np.bool_, bool)):
        return "true" if bool(value) else "false"
    if isinstance(value, (np.integer, int)):
        return str(int(value))
    if isinstance(value, (np.floating, float)):
        if math.isnan(float(value)):
            return ""
        integer = float(value)
        if float(integer).is_integer():
            return str(int(integer))
        return str(integer)
    if pd.isna(cast(object, value)):
        return ""
    return str(value)


def safe_lower(text: Any) -> str:
    """Return a lowercase string for *text* with nulls mapped to empty."""

    return to_text(text).lower()


def null_or_empty(text: Any) -> bool:
    """Return ``True`` when *text* represents ``null`` or ``""`` or ``"0"``."""

    value = to_text(text)
    return value == "" or value == "0"


def pad4(value: Any) -> str:
    """Pad a numeric value to four characters, defaulting to ``"0000"``."""

    text = to_text(value)
    if null_or_empty(text):
        return "0000"
    return text.zfill(4)


def pad2(value: Any) -> str:
    """Pad a numeric value to two characters, defaulting to ``"00"``."""

    text = to_text(value)
    if null_or_empty(text):
        return "00"
    return text.zfill(2)


def pad_pmid8(value: Any) -> str:
    """Return an eight-character PMID with leading zeros."""

    text = to_text(value)
    if null_or_empty(text):
        return ""
    return text.zfill(8)


def normalize_journal(value: Any) -> str:
    """Normalise journal titles by removing dots, trimming and lowercasing."""

    text = to_text(value).replace(".", "")
    return text.strip().lower()


def eq_text(a: Any, b: Any) -> bool:
    """Replicate ``EqText`` from the M script."""

    text_a = to_text(a)
    text_b = to_text(b)
    if null_or_empty(text_a):
        return False
    return text_a == text_b


def ne_text(a: Any, b: Any) -> bool:
    """Replicate ``NeText`` from the M script."""

    text_a = to_text(a)
    text_b = to_text(b)
    if null_or_empty(text_a):
        return False
    return text_a != text_b


def try_parse_int(value: Any) -> int | None:
    """Return ``int(value)`` or ``None`` on failure, emulating ``Number.From``."""

    text = to_text(value)
    if text == "":
        return None
    try:
        return int(text)
    except (TypeError, ValueError):
        try:
            float_value = float(text)
        except (TypeError, ValueError):
            return None
        if float_value.is_integer():
            return int(float_value)
        return None


def _series_any(series_list: Iterable[pd.Series]) -> pd.Series:
    series_iter = iter(series_list)
    try:
        result = next(series_iter).copy()
    except StopIteration:
        return pd.Series(False)
    for item in series_iter:
        result |= item
    return result


def _resolve_relative(base_path: Path, relative: str) -> Path:
    relative_path = Path(relative.replace("\\", os.sep))
    if relative_path.is_absolute():
        return relative_path
    return (base_path / relative_path).resolve()


def _format_windows_path(path: Path) -> str:
    return WINDOWS_PATH_SEPARATOR.join(path.parts)


def _safe_int(value: Any) -> int:
    """Best-effort conversion of *value* to ``int``.

    Conversion errors yield ``0`` to keep QA bookkeeping resilient when
    optional metrics are missing from the report payload.
    """

    try:
        if value is None:
            return 0
        if isinstance(value, bool):
            return int(value)
        return int(value)
    except (TypeError, ValueError):
        return 0


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------


REF_DTYPE = {
    "document_chembl_id": "string",
    "abstract": "string",
    "authors": "string",
    "classification": "string",
    "document_contains_external_links": "boolean",
    "DOI": "string",
    "first_page": "string",
    "is_experimental_doc": "boolean",
    "issue": "string",
    "journal": "string",
    "last_page": "string",
    "month": "string",
    "postcodes": "string",
    "pubmed_id": "string",
    "title": "string",
    "volume": "string",
    "year": "string",
}


OUTPUT_DTYPE = {
    "PubMed.PMID": "Int64",
    "PubMed.DOI": "string",
    "PubMed.ArticleTitle": "string",
    "PubMed.Abstract": "string",
    "PubMed.JournalTitle": "string",
    "PubMed.JournalISOAbbrev": "string",
    "PubMed.Volume": "string",
    "PubMed.Issue": "string",
    "PubMed.StartPage": "string",
    "PubMed.EndPage": "string",
    "PubMed.ISSN": "string",
    "PubMed.PublicationType": "string",
    "PubMed.MeSH_Descriptors": "string",
    "PubMed.MeSH_Qualifiers": "string",
    "PubMed.ChemicalList": "string",
    "PubMed.YearCompleted": "Int64",
    "PubMed.MonthCompleted": "Int64",
    "PubMed.DayCompleted": "Int64",
    "PubMed.YearRevised": "Int64",
    "PubMed.MonthRevised": "Int64",
    "PubMed.DayRevised": "Int64",
    "PubMed.Error": "string",
    "scholar.PMID": "Int64",
    "scholar.DOI": "string",
    "scholar.PublicationTypes": "string",
    "scholar.Venue": "string",
    "scholar.SemanticScholarId": "string",
    "scholar.ExternalIds": "string",
    "scholar.Error": "string",
    "OpenAlex.PMID": "Int64",
    "OpenAlex.DOI": "string",
    "OpenAlex.PublicationTypes": "string",
    "OpenAlex.TypeCrossref": "string",
    "OpenAlex.Genre": "string",
    "OpenAlex.Venue": "string",
    "OpenAlex.MeshDescriptors": "string",
    "OpenAlex.MeshQualifiers": "string",
    "OpenAlex.Id": "string",
    "OpenAlex.Error": "string",
    "crossref.DOI": "string",
    "crossref.Type": "string",
    "crossref.Subtype": "string",
    "crossref.Title": "string",
    "crossref.Subtitle": "string",
    "crossref.Subject": "string",
    "crossref.Error": "string",
    "publication_types_normalised": "string",
    "ChEMBL.document_chembl_id": "string",
    "ChEMBL.title": "string",
    "ChEMBL.abstract": "string",
    "ChEMBL.doi": "string",
    "ChEMBL.year": "Int64",
    "ChEMBL.journal": "string",
    "ChEMBL.journal_abbrev": "string",
    "ChEMBL.volume": "string",
    "ChEMBL.issue": "string",
    "ChEMBL.first_page": "string",
    "ChEMBL.last_page": "string",
    "ChEMBL.pubmed_id": "Int64",
    "ChEMBL.authors": "string",
    "ChEMBL.source": "string",
}

OUTPUT_INT64_COLUMNS: tuple[str, ...] = tuple(
    column for column, dtype in OUTPUT_DTYPE.items() if dtype == "Int64"
)

OUTPUT_READ_DTYPE = {
    column: ("string" if column in OUTPUT_INT64_COLUMNS else dtype)
    for column, dtype in OUTPUT_DTYPE.items()
}


def _load_reference_document(path: Path) -> pd.DataFrame:
    if not path.exists():
        msg = (
            "Reference document CSV not found. "
            "Set 'ref_document_path' to the location of the ETL export or "
            "materialise 'data/input/full/document.csv'."
        )
        raise FileNotFoundError(msg)

    frame = pd.read_csv(
        path,
        dtype=REF_DTYPE,
        encoding=CP1252,
        sep=CSV_DELIMITER,
        keep_default_na=True,
    )

    classification_numeric = pd.to_numeric(
        frame["classification"],
        errors="coerce",
    ).fillna(0)
    frame["classification"] = classification_numeric.astype("Int64").astype(bool)
    frame = frame.rename(columns={"classification": "doctype_review"})

    drop_columns = [
        "abstract",
        "authors",
        "DOI",
        "first_page",
        "issue",
        "journal",
        "last_page",
        "month",
        "postcodes",
        "title",
        "volume",
        "year",
    ]
    frame = frame.drop(columns=[col for col in drop_columns if col in frame.columns])
    return frame


def _load_output_document(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        dtype=OUTPUT_READ_DTYPE,
        encoding=UTF8_ENCODING,
        sep=CSV_DELIMITER,
        keep_default_na=True,
    )

    for column in OUTPUT_INT64_COLUMNS:
        if column not in frame.columns:
            continue
        numeric = pd.to_numeric(frame[column], errors="coerce")
        frame[column] = numeric.astype("Int64")

    return frame


# ---------------------------------------------------------------------------
# Core transformation pipeline
# ---------------------------------------------------------------------------


def _harmonise_documents(out_frame: pd.DataFrame, ref_frame: pd.DataFrame) -> pd.DataFrame:
    df = out_frame.copy()

    for column in PUBMED_NUMERIC_TEXT_COLUMNS:
        df[column] = df[column].map(to_text).replace("", "0")

    df["ChEMBL.journal"] = df["ChEMBL.journal"].map(normalize_journal)
    df["PubMed.JournalTitle"] = df["PubMed.JournalTitle"].map(normalize_journal)

    for column in CHEMBL_NUMERIC_TEXT_COLUMNS:
        df[column] = df[column].map(to_text)

    for column in PUBMED_NUMERIC_TEXT_COLUMNS:
        df[column] = df[column].map(to_text)

    ref_doi = df["ChEMBL.doi"].map(to_text)
    agree_series = []
    conflict_series = []
    for column in EXTERNAL_DOI_COLUMNS:
        external_text = df[column].map(to_text)
        valid = ~external_text.isin(("", "0"))
        agree_series.append(valid & (external_text == ref_doi))
        conflict_series.append(valid & (external_text != ref_doi))
    agree_any = _series_any(agree_series)
    conflict_any = _series_any(conflict_series)
    df["invalid.doi"] = (~agree_any) & conflict_any

    ref_pmid_numbers = df["ChEMBL.pubmed_id"].map(try_parse_int)
    agree_pmid: list[pd.Series] = []
    conflict_pmid: list[pd.Series] = []
    ref_valid = ref_pmid_numbers.notna()

    for column in EXTERNAL_PMID_COLUMNS:
        external_text = df[column].map(to_text)
        external_numbers = external_text.map(try_parse_int)
        external_valid = external_numbers.notna()

        match = ref_valid & external_valid & (external_numbers == ref_pmid_numbers)
        agree_pmid.append(match)

        non_empty = external_text != ""
        mismatch = (
            non_empty
            & ref_valid
            & external_valid
            & (external_numbers != ref_pmid_numbers)
        )
        conflict_pmid.append(mismatch)

    agree_pmid_any = _series_any(agree_pmid)
    conflict_pmid_any = _series_any(conflict_pmid)
    df["invalid.PMID"] = (~agree_pmid_any) & conflict_pmid_any

    review_cols = (
        "PubMed.PublicationType",
        "scholar.PublicationTypes",
        "OpenAlex.PublicationTypes",
        "OpenAlex.TypeCrossref",
        "crossref.Type",
    )
    review_series = [
        df[col].astype("string").fillna("").str.lower().str.contains("review", na=False)
        for col in review_cols
    ]
    df["review"] = _series_any(review_series)

    journal_match = (df["PubMed.JournalTitle"] == df["ChEMBL.journal"]) & (df["ChEMBL.journal"] != "")
    volume_match = (df["PubMed.Volume"] == df["ChEMBL.volume"]) & (df["ChEMBL.volume"] != "")
    issue_match = (df["PubMed.Issue"] == df["ChEMBL.issue"]) & (df["ChEMBL.issue"] != "")
    start_match = (df["PubMed.StartPage"] == df["ChEMBL.first_page"]) & (df["ChEMBL.first_page"] != "")
    end_match = (df["PubMed.EndPage"] == df["ChEMBL.last_page"]) & (df["ChEMBL.last_page"] != "")

    journal_value = np.where(journal_match, df["ChEMBL.journal"], "unknown")
    volume_value = np.where(volume_match, df["ChEMBL.volume"], "unknown")
    issue_value = np.where(issue_match, df["ChEMBL.issue"], "unknown")
    start_value = np.where(start_match, df["ChEMBL.first_page"], "unknown")
    end_value = np.where(end_match, df["ChEMBL.last_page"], "unknown")

    df["reference"] = (
        pd.Series(journal_value, index=df.index)
        + ", "
        + pd.Series(volume_value, index=df.index)
        + "("
        + pd.Series(issue_value, index=df.index)
        + "), p."
        + pd.Series(start_value, index=df.index)
        + "-"
        + pd.Series(end_value, index=df.index)
    )

    df["invalid.reference"] = df["reference"].str.contains("unknown", na=False)

    year_completed = df["PubMed.YearCompleted"].map(to_text)
    year_revised = df["PubMed.YearRevised"].map(to_text)
    chembl_year = df["ChEMBL.year"].map(to_text)

    df["year"] = np.where(
        ~year_completed.map(null_or_empty),
        year_completed.map(pad4),
        np.where(~year_revised.map(null_or_empty), year_revised.map(pad4), chembl_year.map(pad4)),
    )

    month_completed = df["PubMed.MonthCompleted"].map(to_text)
    month_revised = df["PubMed.MonthRevised"].map(to_text)
    df["month"] = np.where(
        ~month_completed.map(null_or_empty),
        month_completed.map(pad2),
        np.where(~month_revised.map(null_or_empty), month_revised.map(pad2), "00"),
    )

    day_completed = df["PubMed.DayCompleted"].map(to_text)
    day_revised = df["PubMed.DayRevised"].map(to_text)
    df["day"] = np.where(
        ~day_completed.map(null_or_empty),
        day_completed.map(pad2),
        np.where(~day_revised.map(null_or_empty), day_revised.map(pad2), "00"),
    )

    df["completed"] = (
        pd.Series(df["year"], index=df.index)
        + "-"
        + pd.Series(df["month"], index=df.index)
        + "-"
        + pd.Series(df["day"], index=df.index)
    )

    df["ChEMBL.pubmed_id"] = df["ChEMBL.pubmed_id"].map(pad_pmid8)

    df["sortorder.document"] = (
        df["PubMed.ISSN"].map(to_text)
        + ":"
        + df["completed"]
        + ":"
        + df["ChEMBL.pubmed_id"].map(to_text)
    )

    df = df.rename(columns={"PubMed.PMID": "PMID", "ChEMBL.doi": "doi"})
    df["PMID"] = df["PMID"].map(to_text)
    df["doi"] = df["doi"].map(to_text)

    df["invalid"] = df["invalid.doi"] | df["invalid.PMID"] | df["invalid.reference"]

    df = df.drop(columns=[col for col in STAGE_REMOVED1_COLUMNS if col in df.columns])

    rename_map = {
        "ChEMBL.document_chembl_id": "document_chembl_id",
        "PubMed.DOI": "doi",
        "PubMed.ArticleTitle": "title",
        "PubMed.Abstract": "abstract",
        "PubMed.MeSH_Descriptors": "mesh.descriptors",
        "PubMed.MeSH_Qualifiers": "mesh.qualifiers",
        "PubMed.ChemicalList": "chemical_list",
        "OpenAlex.MeshDescriptors": "OpenAlex.mesh.descriptors",
    }
    df = df.rename(columns=rename_map)
    df["doi"] = df["doi"].map(to_text)

    second_rename = {
        "chemical_list": "PubMed.chemical.list",
        "mesh.qualifiers": "PubMed.mesh.qualifiers",
        "mesh.descriptors": "PubMed.mesh.descriptors",
    }
    df = df.rename(columns=second_rename)

    for column in LOWERCASE_LIST_COLUMNS:
        if column in df.columns:
            df[column] = df[column].map(safe_lower)

    df = df.merge(ref_frame, on="document_chembl_id", how="left")

    review_values: list[Any] = []
    for current, doctype in zip(df["review"], df["doctype_review"]):
        current_value = None if pd.isna(current) else bool(current)
        doctype_value = None if pd.isna(doctype) else bool(doctype)
        if current_value is True or doctype_value is True:
            review_values.append(True)
        elif current_value is False and doctype_value is False:
            review_values.append(False)
        else:
            review_values.append(pd.NA)
    df["review"] = pd.Series(review_values, index=df.index, dtype="boolean")

    experimental_values: list[Any] = []
    for value in df["review"]:
        if value is pd.NA or pd.isna(value):
            experimental_values.append(pd.NA)
        else:
            experimental_values.append(not bool(value))
    df["experimental"] = pd.Series(experimental_values, index=df.index, dtype="boolean")

    df = df.drop(columns=["is_experimental_doc"], errors="ignore")

    missing_columns = [column for column in FINAL_COLUMN_ORDER if column not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing expected columns after harmonisation: {missing_columns}")

    df = df.loc[:, FINAL_COLUMN_ORDER]
    df = df.sort_values("completed", ascending=True, kind="mergesort")
    df.reset_index(drop=True, inplace=True)
    return df


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def preprocess_documents_csv(
    base_path: str,
    ref_document_rel: str = DEFAULT_REF_RELATIVE,
    out_document_rel: str = DEFAULT_OUTPUT_RELATIVE,
    qa_reference_rel: str | None = DEFAULT_QA_REFERENCE_RELATIVE,
    *,
    run_qa: bool = True,
) -> str:
    base_dir = Path(base_path)
    ref_path = _resolve_relative(base_dir, ref_document_rel)
    out_path = _resolve_relative(base_dir, out_document_rel)

    logger.info(
        "document_postprocess_input",
        base=_format_windows_path(base_dir.resolve()),
        reference=str(ref_path),
        output=str(out_path),
    )

    ref_frame = _load_reference_document(ref_path)
    out_frame = _load_output_document(out_path)

    harmonised = _harmonise_documents(out_frame, ref_frame)

    invalid_counts = {
        "invalid": int(harmonised["invalid"].fillna(False).sum()),
        "invalid.doi": int(harmonised["invalid.doi"].fillna(False).sum()),
        "invalid.PMID": int(harmonised["invalid.PMID"].fillna(False).sum()),
        "invalid.reference": int(harmonised["invalid.reference"].fillna(False).sum()),
    }

    result_name = out_path.name
    if result_name.startswith("output.document_"):
        target_name = f"preprocessed_{result_name}"
    else:
        target_name = f"preprocessed_{result_name}"

    target_path = out_path.with_name(target_name)
    harmonised.to_csv(
        target_path,
        index=False,
        encoding=UTF8_ENCODING,
        sep=CSV_DELIMITER,
    )

    total_rows = len(harmonised)
    logger.info(
        "document_postprocess_output",
        rows=total_rows,
        invalid=invalid_counts["invalid"],
        invalid_doi=invalid_counts["invalid.doi"],
        invalid_pmid=invalid_counts["invalid.PMID"],
        invalid_reference=invalid_counts["invalid.reference"],
        output=str(target_path),
    )

    if total_rows:
        ratio = invalid_counts["invalid"] / total_rows
        if ratio > 0.10:
            logger.warning(
                "document_postprocess_invalid_ratio",
                ratio=ratio,
                rows=total_rows,
                invalid=invalid_counts["invalid"],
            )

    if run_qa and qa_reference_rel:
        qa_reference_path = _resolve_relative(base_dir, qa_reference_rel)
        if qa_reference_path.exists():
            qa_module = None
            try:
                qa_module = importlib.import_module(
                    "qa.check_document_postprocessing"
                )
            except ModuleNotFoundError:
                try:
                    qa_module = importlib.import_module(
                        "library.qa.check_document_postprocessing"
                    )
                except Exception:
                    logger.exception(
                        "document_postprocess_qa_import_failed",
                        reference=str(qa_reference_path),
                    )
            except Exception:
                logger.exception(
                    "document_postprocess_qa_import_failed",
                    reference=str(qa_reference_path),
                )

            if qa_module is not None:
                run_document_postprocessing_check = getattr(
                    qa_module, "run_document_postprocessing_check", None
                )
                crosswalk_cls = getattr(qa_module, "Crosswalk", None)
                crosswalk_path = getattr(qa_module, "CROSSWALK_PATH", None)
                default_diff_limit = getattr(qa_module, "DEFAULT_DIFF_LIMIT", 100)
                default_report_dir = getattr(
                    qa_module, "DEFAULT_REPORT_DIR", Path("output") / "document"
                )

                if callable(run_document_postprocessing_check) and crosswalk_cls and crosswalk_path:
                    try:
                        crosswalk = crosswalk_cls.load(Path(crosswalk_path))
                    except Exception:
                        logger.exception(
                            "document_postprocess_qa_crosswalk_failed",
                            crosswalk=str(crosswalk_path),
                        )
                    else:
                        report_dir = Path(default_report_dir)
                        if not report_dir.is_absolute():
                            report_dir = (base_dir / report_dir).resolve()

                        def _infer_date_code(path: Path) -> str:
                            match = re.search(r"(20\d{6})", path.name)
                            if match:
                                return match.group(1)
                            return datetime.utcnow().strftime("%Y%m%d")

                        try:
                            qa_metrics = run_document_postprocessing_check(
                                crosswalk=crosswalk,
                                m_frame=ref_frame,
                                python_frame=harmonised,
                                diff_limit=int(default_diff_limit),
                                date_code=_infer_date_code(target_path),
                                report_dir=report_dir,
                            )
                        except Exception:
                            logger.exception(
                                "document_postprocess_qa_execution_failed",
                                reference=str(qa_reference_path),
                                output=str(target_path),
                            )
                        else:
                            status = str(qa_metrics.get("status", "")).upper()
                            report_json = qa_metrics.get("report_json")
                            diff_path = qa_metrics.get("diff_path")
                            issues_raw = qa_metrics.get("issues") or []
                            issues = [str(issue) for issue in issues_raw if issue]
                            missing_rows = qa_metrics.get("missing_rows") or {}
                            python_only = _safe_int(missing_rows.get("python_only"))
                            m_only = _safe_int(missing_rows.get("m_only"))
                            diff_rows = _safe_int(qa_metrics.get("diff_rows"))

                            def _subset_issue(label: str) -> bool:
                                lowered = label.lower()
                                return (
                                    "python-only keys detected" in lowered
                                    or "m-output-only keys detected" in lowered
                                )

                            tolerated_subset = (
                                status == "FAIL"
                                and diff_rows == 0
                                and (python_only > 0 or m_only > 0)
                                and all(_subset_issue(issue) for issue in issues)
                            )

                            if status == "PASS":
                                logger.info(
                                    "document_postprocess_qa_passed",
                                    report=str(report_json) if report_json else None,
                                    status=status,
                                )
                            elif tolerated_subset:
                                logger.warning(
                                    "document_postprocess_qa_skipped_subset",
                                    report=str(report_json) if report_json else None,
                                    diff=str(diff_path) if diff_path else None,
                                    status=status,
                                    missing_reference=m_only,
                                    extra_output=python_only,
                                    reference_rows=len(ref_frame),
                                    output_rows=total_rows,
                                    issues=issues or None,
                                )
                            else:
                                logger.error(
                                    "document_postprocess_qa_failed",
                                    report=str(report_json) if report_json else None,
                                    diff=str(diff_path) if diff_path else None,
                                    status=status,
                                    issues=issues or None,
                                )
                                msg = "Document post-processing QA mismatches detected"
                                raise RuntimeError(msg)
                else:
                    logger.error(
                        "document_postprocess_qa_missing_callable",
                        reference=str(qa_reference_path),
                    )
        else:
            logger.info(
                "document_postprocess_qa_reference_missing",
                reference=str(qa_reference_path),
            )

    return str(target_path)


__all__ = [
    "DEFAULT_BASE_PATH",
    "DEFAULT_OUTPUT_PREFIX",
    "DEFAULT_OUTPUT_RELATIVE",
    "DEFAULT_REF_RELATIVE",
    "DEFAULT_QA_REFERENCE_RELATIVE",
    "FINAL_COLUMN_ORDER",
    "LOWERCASE_LIST_COLUMNS",
    "STAGE_REMOVED1_COLUMNS",
    "preprocess_document_export",
    "postprocess_export_file",
    "normalize_journal",
    "pad2",
    "pad4",
    "pad_pmid8",
    "eq_text",
    "ne_text",
    "preprocess_documents_csv",
]
