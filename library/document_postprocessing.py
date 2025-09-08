"""Post-processing utilities for document metadata.

This module normalises and enriches the combined document information
produced by :mod:`get_document_data.py`.  The transformation is a
translation of a Power Query script into ``pandas`` operations.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

import pandas as pd

logger = logging.getLogger(__name__)

# Columns that should be treated as text
TEXT_COLUMNS: list[str] = [
    "document_chembl_id",
    "title",
    "abstract",
    "doi",
    "journal",
    "journal_abbrev",
    "authors",
    "source",
    "PubMed.DOI",
    "PubMed.ArticleTitle",
    "PubMed.Abstract",
    "PubMed.JournalTitle",
    "PubMed.ISSN",
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
    "scholar.Venue",
    "scholar.PublicationTypes",
    "scholar.SemanticScholarId",
    "scholar.ExternalIds",
    "scholar.DOI",
    "scholar.Error",
    "OpenAlex.PublicationTypes",
    "OpenAlex.TypeCrossref",
    "OpenAlex.Genre",
    "OpenAlex.Id",
    "OpenAlex.Venue",
    "OpenAlex.MeshDescriptors",
    "OpenAlex.MeshQualifiers",
    "OpenAlex.Error",
    "crossref.Type",
    "crossref.Subtype",
    "crossref.Title",
    "crossref.Subtitle",
    "crossref.Subject",
    "crossref.Error",
]

# Columns that should be converted to integers
INT_COLUMNS: list[str] = [
    "year",
    "volume",
    "issue",
    "first_page",
    "last_page",
    "pubmed_id",
    "PubMed.PMID",
    "PubMed.Volume",
    "PubMed.Issue",
    "PubMed.StartPage",
    "PubMed.EndPage",
    "scholar.PMID",
]

# Final column ordering of the exported table
COLUMN_ORDER: list[str] = [
    "document_chembl_id",
    "title",
    "abstract",
    "doi",
    "year",
    "journal",
    "journal_abbrev",
    "volume",
    "issue",
    "first_page",
    "last_page",
    "pubmed_id",
    "authors",
    "source",
    "PubMed.PMID",
    "PubMed.DOI",
    "PubMed.ArticleTitle",
    "PubMed.Abstract",
    "PubMed.JournalTitle",
    "PubMed.Volume",
    "PubMed.Issue",
    "PubMed.StartPage",
    "PubMed.EndPage",
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
    "scholar.PMID",
    "scholar.Venue",
    "scholar.SemanticScholarId",
    "scholar.ExternalIds",
    "scholar.DOI",
    "scholar.Error",
    "OpenAlex.TypeCrossref",
    "OpenAlex.Genre",
    "OpenAlex.Id",
    "OpenAlex.Venue",
    "OpenAlex.MeshDescriptors",
    "OpenAlex.MeshQualifiers",
    "OpenAlex.Error",
    "crossref.Type",
    "crossref.Subtype",
    "crossref.Title",
    "crossref.Subtitle",
    "crossref.Subject",
    "crossref.Error",
    "PubMed.ISSN",
    "date_code",
    "Index",
    "PubMed.is_review",
    "scholar.is_review",
    "OpenAlex.is_review",
]


DATE_COLUMNS: Iterable[str] = [
    "PubMed.DayCompleted",
    "PubMed.MonthCompleted",
    "PubMed.YearCompleted",
    "PubMed.DayRevised",
    "PubMed.MonthRevised",
    "PubMed.YearRevised",
]


def _safe_numeric(series: pd.Series) -> pd.Series:
    """Return ``series`` as ``Int64`` with invalid values set to ``<NA>``."""

    return pd.to_numeric(series, errors="coerce").astype("Int64")


def postprocess_documents(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and enrich document metadata.

    Parameters
    ----------
    df:
        Combined table produced by :mod:`get_document_data.py`.

    Returns
    -------
    pandas.DataFrame
        Normalised table with derived fields.
    """

    result = df.copy()
    if result.empty:
        return result

    # --- type conversions -----------------------------------------------------
    for col in TEXT_COLUMNS:
        if col in result.columns:
            result[col] = result[col].astype(str)
    for col in INT_COLUMNS:
        if col in result.columns:
            result[col] = _safe_numeric(result[col])

    for col in DATE_COLUMNS:
        if col not in result.columns:
            result[col] = ""
        else:
            result[col] = result[col].astype(str)

    # --- compute date_code ----------------------------------------------------
    def build_date_code(row: pd.Series) -> str:
        day_completed = row["PubMed.DayCompleted"]
        if day_completed and day_completed != "0":
            return f"{row['PubMed.YearCompleted']}-{row['PubMed.MonthCompleted']}-{day_completed}"
        day_revised = row["PubMed.DayRevised"]
        if day_revised and day_revised != "0":
            return f"{row['PubMed.YearRevised']}-{row['PubMed.MonthRevised']}-{day_revised}"
        return f"{row['PubMed.YearCompleted']}-{row['PubMed.MonthCompleted']}-01"

    result["date_code"] = result.apply(build_date_code, axis=1)

    # --- sort and add running index ------------------------------------------
    result = result.sort_values(
        "date_code", ascending=True, na_position="last"
    ).reset_index(drop=True)
    result["Index"] = result.index

    # --- review flags ---------------------------------------------------------
    def contains_review(series: pd.Series) -> pd.Series:
        return series.str.contains("review", case=False, na=False)

    if "OpenAlex.PublicationTypes" in result.columns:
        result["OpenAlex.is_review"] = contains_review(
            result["OpenAlex.PublicationTypes"]
        )
        result = result.drop(columns=["OpenAlex.PublicationTypes"])
    else:
        result["OpenAlex.is_review"] = False

    if "PubMed.PublicationType" in result.columns:
        result["PubMed.is_review"] = contains_review(result["PubMed.PublicationType"])
        result = result.drop(columns=["PubMed.PublicationType"])
    else:
        result["PubMed.is_review"] = False

    if "scholar.PublicationTypes" in result.columns:
        result["scholar.is_review"] = contains_review(
            result["scholar.PublicationTypes"]
        )
        result = result.drop(columns=["scholar.PublicationTypes"])
    else:
        result["scholar.is_review"] = False

    # --- final column ordering ------------------------------------------------
    existing = [c for c in COLUMN_ORDER if c in result.columns]
    result = result[existing]

    # --- format index ---------------------------------------------------------
    if "Index" in result.columns:
        result["Index"] = result["Index"].astype(int).astype(str).str.zfill(4)

    return result


def postprocess_file(
    input_path: Path | str,
    output_path: Path | str,
    *,
    sep: str = ",",
    encoding: str = "utf8",
) -> None:
    """Read a CSV, apply :func:`postprocess_documents` and write result.

    Parameters
    ----------
    input_path:
        CSV file produced by ``get_document_data.py all``.
    output_path:
        Destination for the cleaned CSV file.
    sep:
        Field delimiter of the CSV files.
    encoding:
        Text encoding of the CSV files.
    """

    df = pd.read_csv(input_path, sep=sep, encoding=encoding, dtype=str)
    processed = postprocess_documents(df)
    processed.to_csv(output_path, index=False, sep=sep, encoding=encoding)
