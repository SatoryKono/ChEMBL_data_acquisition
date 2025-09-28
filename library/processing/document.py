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
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from ..chembl_document import DOCUMENT_COLUMNS as _CHEMBL_COLUMNS
from ..document_type_classifier import compute_scores, decide_label
from ..document_type_terms import parse_terms

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

    if isinstance(data, pd.DataFrame):
        df = data.copy()
    else:
        if not data:
            return pd.DataFrame(columns=list(columns))
        df = pd.DataFrame.from_records(data)

    if fill_missing:
        for col in columns:
            if col not in df.columns:
                df[col] = ""

    extra = sorted(c for c in df.columns if c not in columns)
    head = [c for c in columns if c in df.columns]
    ordered = head + extra if not fill_missing else list(columns) + extra
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
