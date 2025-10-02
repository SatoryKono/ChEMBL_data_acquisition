"""Post-processing stage for harmonising document exports.

The implementation mirrors the behaviour of the Power Query (M) script used in
the original ETL pipeline.  The module accepts the exported CSV from
``scripts.get_document_data`` and the reference ChEMBL document catalogue,
applies the exact transformations described in the M script and writes a
deterministic CSV containing the harmonised view.
"""

from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from library.log import logger


# ---------------------------------------------------------------------------
# Parameters (paths, encodings, delimiters)
# ---------------------------------------------------------------------------

CP1252 = "cp1252"
UTF8 = "utf-8"
CSV_DELIMITER = ","

DEFAULT_BASE_PATH = "e:\\github\\ChEMBL_data_acquisition\\data\\"
DEFAULT_REF_RELATIVE = "input\\full\\document.csv"
DEFAULT_OUTPUT_RELATIVE = "output\\document\\output.document_YYYYMMDD.csv"


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

LIST_SEPARATOR = "\\"


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
    if pd.isna(value):  # type: ignore[arg-type]
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
    return LIST_SEPARATOR.join(path.parts)


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------


REF_DTYPE = {
    "document_chembl_id": "string",
    "abstract": "string",
    "authors": "string",
    "classification": "Int64",
    "document_contains_external_links": "boolean",
    "DOI": "string",
    "first_page": "Int64",
    "is_experimental_doc": "boolean",
    "issue": "Int64",
    "journal": "string",
    "last_page": "Int64",
    "month": "Int64",
    "postcodes": "string",
    "pubmed_id": "Int64",
    "title": "string",
    "volume": "Int64",
    "year": "Int64",
}


OUTPUT_DTYPE = {
    "PubMed.PMID": "Int64",
    "PubMed.DOI": "string",
    "PubMed.ArticleTitle": "string",
    "PubMed.Abstract": "string",
    "PubMed.JournalTitle": "string",
    "PubMed.JournalISOAbbrev": "string",
    "PubMed.Volume": "Int64",
    "PubMed.Issue": "Int64",
    "PubMed.StartPage": "Int64",
    "PubMed.EndPage": "Int64",
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
    "ChEMBL.volume": "Int64",
    "ChEMBL.issue": "Int64",
    "ChEMBL.first_page": "Int64",
    "ChEMBL.last_page": "Int64",
    "ChEMBL.pubmed_id": "Int64",
    "ChEMBL.authors": "string",
    "ChEMBL.source": "string",
}


def _load_reference_document(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        dtype=REF_DTYPE,
        encoding=CP1252,
        sep=CSV_DELIMITER,
        keep_default_na=True,
    )

    frame["classification"] = frame["classification"].fillna(0)
    frame["classification"] = frame["classification"].astype("Int64")
    frame["classification"] = frame["classification"].astype(bool)
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
    return pd.read_csv(
        path,
        dtype=OUTPUT_DTYPE,
        encoding=UTF8,
        sep=CSV_DELIMITER,
        keep_default_na=True,
    )


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
    harmonised.to_csv(target_path, index=False, encoding=UTF8, sep=CSV_DELIMITER)

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

    return str(target_path)


__all__ = [
    "DEFAULT_BASE_PATH",
    "DEFAULT_OUTPUT_RELATIVE",
    "DEFAULT_REF_RELATIVE",
    "FINAL_COLUMN_ORDER",
    "LOWERCASE_LIST_COLUMNS",
    "STAGE_REMOVED1_COLUMNS",
    "normalize_journal",
    "pad2",
    "pad4",
    "pad_pmid8",
    "eq_text",
    "ne_text",
    "preprocess_documents_csv",
]
