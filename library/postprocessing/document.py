"""Post-processing helpers for document exports.

Changelog
~~~~~~~~
- 2024-04-??: Initial Python port of the Power Query pipeline used to derive
  analytics-friendly document summaries.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Mapping

import pandas as pd

from ..config import IoCfg
from ..csv_utils import write_csv_deterministic
from ..validation import validate_columns

# ===== Parameters ===========================================================
UTF8_ENCODING = "utf-8"
DEFAULT_OUTPUT_PREFIX = "preprocessed_"
LIST_SEPARATOR = "; "
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
DEFAULT_KEY_COLUMNS: tuple[str, ...] = ("document_chembl_id",)


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
        collected.append(LIST_SEPARATOR.join(terms))
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
        error_tokens.append(LIST_SEPARATOR.join(ordered))
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
    metadata_strings = [LIST_SEPARATOR.join(tokens) for tokens in token_lists]
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
def preprocess_document_export(df: pd.DataFrame) -> pd.DataFrame:
    """Return analytics-oriented projection of ``df``."""

    validate_columns(df, REQUIRED_INPUT_COLUMNS)
    if df.empty:
        result = pd.DataFrame(columns=list(PREPROCESSED_COLUMN_ORDER))
        return result

    result = pd.DataFrame(index=df.index)
    result["document_chembl_id"] = _string_series(df["document_chembl_id"])

    for target, columns in COALESCE_PRIORITIES.items():
        result[target] = _coalesce_columns(df, columns)

    if "publication_class" in df.columns:
        result["publication_class"] = _string_series(df["publication_class"])
    else:
        result["publication_class"] = pd.Series([""] * len(df), index=df.index, dtype="string")

    review_series = pd.Series([False] * len(df), index=df.index)
    if "publication_class" in result.columns:
        review_series = review_series | (
            result["publication_class"].str.lower() == "review"
        )
    for column in REVIEW_COLUMNS:
        if column in df.columns:
            review_series = review_series | _series_to_bool(df[column])
    result["is_review"] = review_series

    metadata_series, metadata_flags = _build_metadata_flags(df)
    result["metadata_sources"] = metadata_series
    result["metadata_source_count"] = metadata_flags.pop("metadata_source_count")
    for column, mask in metadata_flags.items():
        result[column] = mask

    error_sources, has_error = _build_error_sources(df)
    result["error_sources"] = error_sources
    result["has_error"] = has_error

    result["coverage_score"] = result["metadata_source_count"].astype("Int64")
    status_values = [
        _coverage_status(int(score) if pd.notna(score) else 0, bool(error))
        for score, error in zip(
            result["coverage_score"].fillna(0), result["has_error"].fillna(False)
        )
    ]
    result["coverage_status"] = pd.Series(status_values, index=df.index, dtype="string")

    mesh_columns = [column for column in MESH_TERM_COLUMNS if column in df.columns]
    result["mesh_terms"] = _aggregate_terms(df, mesh_columns)

    if "publication_year" in result.columns:
        publication_year = result["publication_year"].replace("", pd.NA)
        result["publication_year"] = pd.to_numeric(
            publication_year, errors="coerce"
        ).astype("Int64")
    else:
        result["publication_year"] = pd.Series([pd.NA] * len(df), index=df.index, dtype="Int64")

    for column in PASS_THROUGH_COLUMNS:
        if column not in df.columns:
            if column == "Index":
                result[column] = pd.Series([""] * len(df), index=df.index, dtype="string")
            else:
                result[column] = pd.Series([""] * len(df), index=df.index, dtype="string")
            continue
        if column == "Index":
            index_series = df[column]
            if pd.api.types.is_numeric_dtype(index_series):
                padded = index_series.fillna(0).astype(int).astype(str)
            else:
                padded = _string_series(index_series)
            result[column] = padded.str.zfill(INDEX_PAD_WIDTH)
        else:
            result[column] = _string_series(df[column])

    existing_columns = [
        column for column in PREPROCESSED_COLUMN_ORDER if column in result.columns
    ]
    return result.loc[:, existing_columns]


def postprocess_export_file(
    input_path: Path | str,
    *,
    cfg: IoCfg,
    output_path: Path | None = None,
    sep: str | None = None,
    encoding: str | None = None,
) -> Path:
    """Read ``input_path``, project the analytics table and write to disk."""

    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding or UTF8_ENCODING
    frame = pd.read_csv(
        input_path,
        sep=sep,
        encoding=encoding,
        dtype=str,
        keep_default_na=False,
    )
    processed = preprocess_document_export(frame)
    destination = Path(output_path) if output_path else Path(input_path).with_name(
        f"{DEFAULT_OUTPUT_PREFIX}{Path(input_path).name}"
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    col_order = [col for col in PREPROCESSED_COLUMN_ORDER if col in processed.columns]
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


# ===== Exports ==============================================================
__all__ = [
    "preprocess_document_export",
    "postprocess_export_file",
]
