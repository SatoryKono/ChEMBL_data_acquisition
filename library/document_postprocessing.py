"""Post-processing utilities for document metadata.

This module normalises and enriches the combined document information
produced by :mod:`scripts.get_document_data`.  The transformation is a
translation of a Power Query script into ``pandas`` operations.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pandas as pd

from . import validation
from .config import IoCfg

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

REQUIRED_COLUMNS = ["PubMed.PMID"]


def _safe_numeric(series: pd.Series) -> pd.Series:
    """Return ``series`` as ``Int64`` with invalid values set to ``<NA>``."""
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def postprocess_documents(
    df: pd.DataFrame, *, required_columns: Iterable[str] | None = None
) -> pd.DataFrame:
    """Clean and enrich document metadata.

    Parameters
    ----------
    df:
        Combined table produced by :mod:`scripts.get_document_data`.
    required_columns:
        Optional columns that must exist in ``df`` before processing. If
        ``None`` (the default), no schema validation is performed.

    Returns
    -------
    pandas.DataFrame
        Normalised table with derived fields.

    """
    if required_columns is not None:
        validation.validate_columns(df, required_columns)
    # Create a view to avoid an expensive deep copy. Subsequent operations
    # either replace columns or return new DataFrames, so the original
    # ``df`` remains unchanged while minimising memory usage.
    result = df.copy(deep=False)
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
    def _normalise_component(value: object, width: int) -> tuple[str, str]:
        """Return raw and padded string representations of a date component."""

        if pd.isna(value):
            raw = ""
        else:
            raw = str(value).strip()
        if raw.isdigit():
            padded = raw.zfill(width)
        else:
            padded = raw
        return raw, padded

    def build_date_code(row: pd.Series) -> str:
        (
            day_completed_raw,
            day_completed,
        ) = _normalise_component(row["PubMed.DayCompleted"], 2)
        _, month_completed = _normalise_component(row["PubMed.MonthCompleted"], 2)
        _, year_completed = _normalise_component(row["PubMed.YearCompleted"], 4)
        day_revised_raw, day_revised = _normalise_component(
            row["PubMed.DayRevised"], 2
        )
        _, month_revised = _normalise_component(row["PubMed.MonthRevised"], 2)
        _, year_revised = _normalise_component(row["PubMed.YearRevised"], 4)

        if day_completed_raw and day_completed_raw != "0":
            return f"{year_completed}-{month_completed}-{day_completed}"
        if day_revised_raw and day_revised_raw != "0":
            return f"{year_revised}-{month_revised}-{day_revised}"
        return f"{year_completed}-{month_completed}-01"

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
    cfg: IoCfg,
    sep: str | None = None,
    encoding: str | None = None,
) -> None:
    """Read a CSV, apply :func:`postprocess_documents` and write result.

    Parameters
    ----------
    input_path:
        CSV file produced by ``scripts/get_document_data.py all``.
    output_path:
        Destination for the cleaned CSV file.
    cfg:
        I/O configuration providing default CSV parameters.
    sep:
        Field delimiter of the CSV files. Defaults to ``cfg.csv_sep``.
    encoding:
        Text encoding of the CSV files. Defaults to ``cfg.csv_encoding``.

    """
    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding
    df = pd.read_csv(input_path, sep=sep, encoding=encoding, dtype=str)
    processed = postprocess_documents(df)
    processed.to_csv(output_path, index=False, sep=sep, encoding=encoding)
