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

from collections import Counter
from collections.abc import Iterable, Iterator, Mapping, Sequence

from pathlib import Path
from typing import Any

import pandas as pd

from .type_classifier import compute_scores, decide_label
from .type_terms import parse_terms
from ...common.pandas_utils import merge_series_prefer_left
from ...schemas.document_spec import (
    DOCUMENT_COLUMN_GROUPS,
    DOCUMENT_SCHEMA_COLUMNS as _DECLARED_SCHEMA_COLUMNS,
)

# ---------------------------------------------------------------------------
CH_EMBL_COLUMNS: list[str] = list(DOCUMENT_COLUMN_GROUPS["chembl"])
DERIVED_COLUMNS: list[str] = list(DOCUMENT_COLUMN_GROUPS["derived"])
PIPELINE_STATUS_COLUMNS: list[str] = list(DOCUMENT_COLUMN_GROUPS["pipeline_status"])
PUBMED_COLUMNS: list[str] = list(DOCUMENT_COLUMN_GROUPS["pubmed"])
SEMANTIC_SCHOLAR_COLUMNS: list[str] = list(DOCUMENT_COLUMN_GROUPS["scholar"])
OPENALEX_COLUMNS: list[str] = list(DOCUMENT_COLUMN_GROUPS["openalex"])
CROSSREF_COLUMNS: list[str] = list(DOCUMENT_COLUMN_GROUPS["crossref"])
RUNTIME_COLUMNS: list[str] = list(DOCUMENT_COLUMN_GROUPS["pipeline_runtime"])

DOCUMENT_SCHEMA_COLUMNS: list[str] = list(_DECLARED_SCHEMA_COLUMNS)

_EXPECTED_ORDER = (
    CH_EMBL_COLUMNS
    + DERIVED_COLUMNS
    + PIPELINE_STATUS_COLUMNS
    + PUBMED_COLUMNS
    + SEMANTIC_SCHOLAR_COLUMNS
    + OPENALEX_COLUMNS
    + CROSSREF_COLUMNS
    + RUNTIME_COLUMNS
)

if _EXPECTED_ORDER != DOCUMENT_SCHEMA_COLUMNS:  # pragma: no cover - config error
    raise RuntimeError(
        "document schema declaration order mismatch between groups and schema"
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
    for prefix in ("doi:", "https://doi.org/", "http://doi.org/", "doi.org/"):
        if doi.lower().startswith(prefix):
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
    chembl_df: pd.DataFrame | Iterable[pd.DataFrame],
    metadata_df: pd.DataFrame | Iterable[pd.DataFrame],
) -> pd.DataFrame:
    """Merge PubMed style metadata into ``chembl_df`` chunk by chunk."""

    def _iter_frames(frame_or_iterable: Iterable[pd.DataFrame]) -> Iterator[pd.DataFrame]:
        return iter(frame_or_iterable)

    if isinstance(chembl_df, pd.DataFrame):
        result = chembl_df.copy()
    else:
        chembl_iter = _iter_frames(chembl_df)
        try:
            result = next(chembl_iter).copy()
        except StopIteration:
            return pd.DataFrame(columns=DOCUMENT_SCHEMA_COLUMNS)
        for extra_frame in chembl_iter:
            if extra_frame.empty:
                continue
            result = pd.concat([result, extra_frame], ignore_index=True)

    if result.empty:
        return result

    join_key = "pubmed_id" if "pubmed_id" in result.columns else "PubMed.PMID"
    if join_key not in result.columns:
        return result

    result[join_key] = (
        pd.to_numeric(result[join_key], errors="coerce")
        .astype("Int64")
        .astype("string")
        .fillna("")
    )

    left_columns = list(result.columns)
    result["__merge_order"] = range(len(result))
    result = result.set_index(join_key, drop=False)

    if isinstance(metadata_df, pd.DataFrame):
        metadata_iter = iter((metadata_df,))
    else:
        metadata_iter = _iter_frames(metadata_df)
    metadata_columns: list[str] = []

    for chunk in metadata_iter:
        if chunk.empty or "PubMed.PMID" not in chunk.columns:
            continue
        right = chunk.drop(columns=[col for col in CH_EMBL_COLUMNS if col in chunk.columns])
        right["PubMed.PMID"] = (
            pd.to_numeric(right["PubMed.PMID"], errors="coerce")
            .astype("Int64")
            .astype("string")
            .fillna("")
        )
        right = right.set_index("PubMed.PMID", drop=False)

        for column in right.columns:
            if column == join_key:
                continue
            if column not in result.columns:
                result[column] = pd.Series(
                    [pd.NA] * len(result), index=result.index, dtype="object"
                )
            if column not in metadata_columns and column not in left_columns:
                metadata_columns.append(column)
            result[column] = merge_series_prefer_left(result[column], right[column])

    result = result.sort_values("__merge_order").drop(columns=["__merge_order"])
    result = result.reset_index(drop=True)

    base_columns = [col for col in left_columns if col in result.columns]
    new_columns = [col for col in metadata_columns if col in result.columns]
    extra_columns = [
        col
        for col in result.columns
        if col not in {"__merge_order", *base_columns, *new_columns}
    ]
    ordered = base_columns + new_columns + extra_columns
    return result[ordered]


def dataframe_to_strings(
    df: pd.DataFrame, *, skip: Iterable[str] | None = None
) -> pd.DataFrame:
    """Convert non-numeric columns of ``df`` to strings."""

    skip_set = set(skip or [])
    result = df.copy()
    dtypes = result.dtypes
    for idx, col in enumerate(result.columns):
        if col in skip_set:
            continue
        if pd.api.types.is_numeric_dtype(dtypes.iloc[idx]):
            continue
        result.iloc[:, idx] = result.iloc[:, idx].astype(str)
    return result


# ---------------------------------------------------------------------------
# Quality reporting


class DocumentQualityAccumulator:
    """Incrementally compute summary metrics for document exports."""

    def __init__(self) -> None:
        self.rows_total = 0
        self._doi_truthy = 0
        self._doi_total = 0
        self._class_counts: Counter[str] = Counter()
        self._error_counts: Counter[str] = Counter()
        self._placeholder_counts: Counter[str] = Counter()

    def consume(self, frame: pd.DataFrame) -> None:
        if not isinstance(frame, pd.DataFrame):
            raise TypeError("DocumentQualityAccumulator.consume expects a DataFrame")

        row_count = len(frame)
        self.rows_total += row_count
        if row_count == 0:
            return

        if "doi" in frame.columns:
            doi_series = frame["doi"].astype(str).str.strip()
        elif "doi_normalised" in frame.columns:
            doi_series = frame["doi_normalised"].astype(str).str.strip()
        else:
            doi_series = None
        if doi_series is not None:
            self._doi_truthy += int((doi_series != "").sum())
            self._doi_total += len(doi_series)

        publication_class = frame.get("publication_class")
        if publication_class is not None:
            cleaned = (
                publication_class.fillna("unknown")
                .astype(str)
                .str.strip()
                .replace("", "unknown")
            )
            self._class_counts.update(cleaned.value_counts().to_dict())

        for column, key in (
            ("PubMed.Error", "pubmed"),
            ("scholar.Error", "semantic_scholar"),
            ("OpenAlex.Error", "openalex"),
            ("crossref.Error", "crossref"),
        ):
            series = frame.get(column)
            if series is None:
                continue
            truthy = (
                series.astype(str)
                .str.strip()
                .astype(bool)
            )
            self._error_counts[key] += int(truthy.sum())

        status = frame.get("fetch_status")
        if status is not None:
            status_strings = status.astype("string").fillna("")
            is_error = status_strings.str.strip().str.lower() == "error"
            if bool(is_error.any()):
                sources = frame.get("error_source")
                if sources is None:
                    source_strings = pd.Series(
                        ["unknown"] * len(frame),
                        index=status_strings.index,
                        dtype="string",
                    )
                else:
                    source_strings = (
                        sources.astype("string")
                        .fillna("unknown")
                        .str.strip()
                        .str.lower()
                    )
                    source_strings = source_strings.replace("", "unknown")
                counts = source_strings[is_error].value_counts()
                self._placeholder_counts.update(
                    {str(source): int(count) for source, count in counts.items()}
                )

    def build(self) -> dict[str, Any]:
        doi_coverage = (
            float(self._doi_truthy / self._doi_total) if self._doi_total else 0.0
        )
        return {
            "rows_total": int(self.rows_total),
            "doi_coverage": doi_coverage,
            "publication_class_counts": dict(self._class_counts),
            "error_counts": {
                "pubmed": int(self._error_counts.get("pubmed", 0)),
                "semantic_scholar": int(
                    self._error_counts.get("semantic_scholar", 0)
                ),
                "openalex": int(self._error_counts.get("openalex", 0)),
                "crossref": int(self._error_counts.get("crossref", 0)),
            },
            "error_placeholder_counts": self._build_placeholder_counts(),
        }

    def _build_placeholder_counts(self) -> dict[str, int]:
        """Return normalised placeholder counters keyed by error source."""

        baseline = {
            "pubmed": 0,
            "semantic_scholar": 0,
            "openalex": 0,
            "crossref": 0,
            "unknown": 0,
        }
        counts = {key: int(self._placeholder_counts.get(key, 0)) for key in baseline}
        for key, value in self._placeholder_counts.items():
            if key not in counts:
                counts[key] = int(value)
        return counts


def build_quality_report(
    df: pd.DataFrame | Iterable[pd.DataFrame] | DocumentQualityAccumulator,
) -> dict[str, Any]:
    """Return a JSON serialisable quality summary for ``df``."""

    if isinstance(df, DocumentQualityAccumulator):
        accumulator = df
    else:
        accumulator = DocumentQualityAccumulator()
        frames: Iterable[pd.DataFrame]
        if isinstance(df, pd.DataFrame):
            frames = (df,)
        elif isinstance(df, Iterable):
            frames = df
        else:
            raise TypeError("Unsupported input for build_quality_report")
        for frame in frames:
            accumulator.consume(frame)

    return accumulator.build()


def save_quality_report(report: Mapping[str, Any], path: Path) -> Path:
    """Persist ``report`` as JSON file and return ``path``."""

    path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return path


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
]
