"""Post-processing utilities for document metadata.

This module implements the Stage_FinalOrder Power Query transformation used
within the legacy ETL pipeline.  The operations are a direct port of the
original M-script and therefore favour explicit helper functions mirroring the
source implementation (``NullOrEmpty``, ``NormalizeJournal`` and friends).
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from ...common.csv_utils import write_csv_deterministic
from ...common.text_utils import to_text
from ...config import IoCfg

# ===== Parameters ===========================================================

UTF8_ENCODING = "utf-8"
CP1252_ENCODING = "cp1252"
CSV_DELIMITER = ","
OUTPUT_PREFIX = "preprocessed_"
DEFAULT_REF_DOCUMENT_PATH = Path("data/input/full/document.csv")


# ===== Column specifications ===============================================

PUBMED_NUMERIC_TEXT_COLUMNS: tuple[str, ...] = (
    "PubMed.Volume",
    "PubMed.Issue",
    "PubMed.StartPage",
    "PubMed.EndPage",
)

CHEMBL_NUMERIC_TEXT_COLUMNS: tuple[str, ...] = (
    "ChEMBL.volume",
    "ChEMBL.issue",
    "ChEMBL.first_page",
    "ChEMBL.last_page",
)

EXTERNAL_DOI_COLUMNS: tuple[str, ...] = (
    "PubMed.DOI",
    "scholar.DOI",
    "OpenAlex.DOI",
    "crossref.DOI",
)

EXTERNAL_PMID_COLUMNS: tuple[str, ...] = (
    "PubMed.PMID",
    "scholar.PMID",
    "OpenAlex.PMID",
)

LOWERCASE_LIST_COLUMNS: tuple[str, ...] = (
    "PubMed.mesh.descriptors",
    "PubMed.mesh.qualifiers",
    "PubMed.chemical.list",
    "OpenAlex.mesh.descriptors",
)

STAGE_REMOVED_COLUMNS: tuple[str, ...] = (
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

FINAL_COLUMN_ORDER: tuple[str, ...] = (
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

REFERENCE_REQUIRED_COLUMNS: tuple[str, ...] = (
    "document_chembl_id",
    "doctype_review",
    "document_contains_external_links",
    "is_experimental_doc",
)


# ===== Helper utilities =====================================================


def safe_lower(value: Any) -> str:
    """Return a lowercase representation with nulls mapped to ``""``."""

    return to_text(value).lower()


def null_or_empty(value: Any) -> bool:
    """Return ``True`` when *value* is ``null``, ``""`` or the literal ``"0"``."""

    text = to_text(value)
    return text in {"", "0"}


def pad2(value: Any) -> str:
    """Pad numbers to two characters using leading zeros."""

    text = to_text(value)
    if null_or_empty(text):
        return "00"
    return text.zfill(2)


def pad4(value: Any) -> str:
    """Pad numbers to four characters using leading zeros."""

    text = to_text(value)
    if null_or_empty(text):
        return "0000"
    return text.zfill(4)


def pad_pmid8(value: Any) -> str:
    """Return an eight-character PMID with leading zeros."""

    text = to_text(value)
    if null_or_empty(text):
        return ""
    return text.zfill(8)


def normalize_journal(value: Any) -> str:
    """Normalise journal titles by removing dots, trimming and lowercasing."""

    return to_text(value).replace(".", "").strip().lower()


def eq_text(a: Any, b: Any) -> bool:
    """Replicate the Power Query ``EqText`` guard."""

    text_a = to_text(a)
    if text_a == "":
        return False
    return text_a == to_text(b)


def ne_text(a: Any, b: Any) -> bool:
    """Replicate the Power Query ``NeText`` guard."""

    text_a = to_text(a)
    if text_a == "":
        return False
    return text_a != to_text(b)


def try_parse_int(value: Any) -> int | None:
    """Return an integer representation or ``None`` on failure."""

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
    iterator = iter(series_list)
    try:
        result = next(iterator).copy()
    except StopIteration:
        return pd.Series(False)
    for series in iterator:
        result |= series
    return result


def _coerce_to_bool(value: Any) -> Any:
    if pd.isna(cast(object, value)):
        return pd.NA
    text = to_text(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    return pd.NA


def _ensure_text_column(df: pd.DataFrame, column: str, default: str = "") -> pd.Series:
    if column in df.columns:
        series = df[column].map(to_text)
    else:
        series = pd.Series(default, index=df.index, dtype="string")
    df[column] = series
    return series


# ===== Reference document handling =========================================


def _normalise_reference_frame(frame: pd.DataFrame) -> pd.DataFrame:
    reference = frame.copy()
    reference = reference.rename(columns={"classification": "doctype_review"})

    if "doctype_review" in reference.columns:
        numeric = pd.to_numeric(reference["doctype_review"], errors="coerce")
        reference["doctype_review"] = numeric.fillna(0).astype("Int64").astype(bool)
    else:
        reference["doctype_review"] = False

    bool_columns = ("document_contains_external_links", "is_experimental_doc")
    for column in bool_columns:
        if column in reference.columns:
            coerced = reference[column].map(_coerce_to_bool)
            reference[column] = coerced.fillna(False).astype("boolean")
        else:
            reference[column] = pd.Series(False, index=reference.index, dtype="boolean")

    if "document_chembl_id" not in reference.columns:
        reference["document_chembl_id"] = pd.Series(dtype="string")

    missing = [
        col for col in REFERENCE_REQUIRED_COLUMNS if col not in reference.columns
    ]
    if missing:
        for column in missing:
            if column == "document_chembl_id":
                reference[column] = pd.Series(dtype="string")
            else:
                reference[column] = pd.Series(
                    False, index=reference.index, dtype="boolean"
                )

    return reference.loc[:, REFERENCE_REQUIRED_COLUMNS]


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
        dtype={
            "document_chembl_id": "string",
            "classification": "string",
            "document_contains_external_links": "boolean",
            "is_experimental_doc": "boolean",
        },
        encoding=CP1252_ENCODING,
        sep=CSV_DELIMITER,
        keep_default_na=True,
    )

    return _normalise_reference_frame(frame)


def _resolve_reference(
    ref_document: pd.DataFrame | None,
    ref_document_path: Path | str | None,
) -> pd.DataFrame:
    if ref_document is not None:
        return _normalise_reference_frame(ref_document)
    if ref_document_path is not None:
        return _load_reference_document(Path(ref_document_path))
    return _load_reference_document(DEFAULT_REF_DOCUMENT_PATH)


# ===== Core transformation ==================================================


def _prepare_input_frame(df: pd.DataFrame) -> pd.DataFrame:
    frame = df.copy()

    if (
        "ChEMBL.document_chembl_id" not in frame.columns
        and "document_chembl_id" in frame.columns
    ):
        frame = frame.rename(
            columns={"document_chembl_id": "ChEMBL.document_chembl_id"}
        )

    if "ChEMBL.document_chembl_id" not in frame.columns:
        raise ValueError("Missing required column 'ChEMBL.document_chembl_id'")

    frame["ChEMBL.document_chembl_id"] = frame["ChEMBL.document_chembl_id"].map(to_text)

    for column in PUBMED_NUMERIC_TEXT_COLUMNS:
        series = _ensure_text_column(frame, column)
        frame[column] = series.replace("", "0")

    for column in CHEMBL_NUMERIC_TEXT_COLUMNS:
        _ensure_text_column(frame, column)

    for column in (
        "ChEMBL.doi",
        "ChEMBL.journal",
        "PubMed.JournalTitle",
        "PubMed.ISSN",
        "ChEMBL.year",
        "ChEMBL.pubmed_id",
        "PubMed.YearCompleted",
        "PubMed.MonthCompleted",
        "PubMed.DayCompleted",
        "PubMed.YearRevised",
        "PubMed.MonthRevised",
        "PubMed.DayRevised",
    ):
        _ensure_text_column(frame, column)

    list_columns = (
        "PubMed.MeSH_Descriptors",
        "PubMed.MeSH_Qualifiers",
        "PubMed.ChemicalList",
        "OpenAlex.MeshDescriptors",
    )
    for column in list_columns:
        _ensure_text_column(frame, column)

    for column in EXTERNAL_DOI_COLUMNS + EXTERNAL_PMID_COLUMNS:
        _ensure_text_column(frame, column)

    review_columns = (
        "PubMed.PublicationType",
        "scholar.PublicationTypes",
        "OpenAlex.PublicationTypes",
        "OpenAlex.TypeCrossref",
        "crossref.Type",
    )
    for column in review_columns:
        _ensure_text_column(frame, column)

    optional_text = (
        "PubMed.ArticleTitle",
        "PubMed.Abstract",
        "scholar.Venue",
        "scholar.SemanticScholarId",
        "scholar.ExternalIds",
        "OpenAlex.Venue",
        "OpenAlex.MeshQualifiers",
        "OpenAlex.Id",
        "crossref.Title",
        "ChEMBL.title",
        "ChEMBL.abstract",
        "ChEMBL.journal_abbrev",
        "ChEMBL.authors",
        "ChEMBL.source",
        "PubMed.Error",
    )
    for column in optional_text:
        _ensure_text_column(frame, column)

    return frame


def postprocess_documents(
    df: pd.DataFrame,
    *,
    required_columns: Iterable[str] | None = None,
    ref_document: pd.DataFrame | None = None,
    ref_document_path: Path | str | None = None,
) -> pd.DataFrame:
    """Clean and enrich document metadata.

    Parameters
    ----------
    df:
        Combined table produced by :mod:`scripts.get_document_data`.
    required_columns:
        Optional schema guard applied before processing.
    ref_document:
        Optional pre-loaded reference table (``ref_document`` from the legacy
        pipeline). When omitted, :data:`DEFAULT_REF_DOCUMENT_PATH` is used.
    ref_document_path:
        Path to the reference CSV. Ignored when ``ref_document`` is supplied.
    """

    if required_columns is not None:
        from ... import validation as validation_module

        validation_module.validate_columns(df, required_columns)

    if df.empty:
        return pd.DataFrame(columns=FINAL_COLUMN_ORDER)

    frame = _prepare_input_frame(df)
    reference = _resolve_reference(ref_document, ref_document_path)

    frame["ChEMBL.journal"] = frame["ChEMBL.journal"].map(normalize_journal)
    frame["PubMed.JournalTitle"] = frame["PubMed.JournalTitle"].map(normalize_journal)

    ref_doi = frame["ChEMBL.doi"].map(to_text)
    doi_agree: list[pd.Series] = []
    doi_conflict: list[pd.Series] = []
    for column in EXTERNAL_DOI_COLUMNS:
        external = frame[column].map(to_text)
        valid = ~external.isin(("", "0"))
        doi_agree.append(valid & (external == ref_doi))
        doi_conflict.append(valid & (external != ref_doi))
    frame["invalid.doi"] = (~_series_any(doi_agree)) & _series_any(doi_conflict)

    ref_pmid_numbers = frame["ChEMBL.pubmed_id"].map(try_parse_int)
    pmid_agree: list[pd.Series] = []
    pmid_conflict: list[pd.Series] = []
    ref_valid = ref_pmid_numbers.notna()
    for column in EXTERNAL_PMID_COLUMNS:
        external_text = frame[column].map(to_text)
        external_numbers = external_text.map(try_parse_int)
        external_valid = external_numbers.notna()
        match = ref_valid & external_valid & (external_numbers == ref_pmid_numbers)
        pmid_agree.append(match)

        non_empty = external_text != ""
        mismatch = (
            non_empty
            & ref_valid
            & external_valid
            & (external_numbers != ref_pmid_numbers)
        )
        pmid_conflict.append(mismatch)
    frame["invalid.PMID"] = (~_series_any(pmid_agree)) & _series_any(pmid_conflict)

    review_columns = (
        "PubMed.PublicationType",
        "scholar.PublicationTypes",
        "OpenAlex.PublicationTypes",
        "OpenAlex.TypeCrossref",
        "crossref.Type",
    )
    review_series = [
        frame[column]
        .astype("string")
        .fillna("")
        .str.lower()
        .str.contains("review", na=False)
        for column in review_columns
    ]
    frame["review"] = _series_any(review_series)

    journal_match = (frame["PubMed.JournalTitle"] == frame["ChEMBL.journal"]) & (
        frame["ChEMBL.journal"] != ""
    )
    volume_match = (frame["PubMed.Volume"] == frame["ChEMBL.volume"]) & (
        frame["ChEMBL.volume"] != ""
    )
    issue_match = (frame["PubMed.Issue"] == frame["ChEMBL.issue"]) & (
        frame["ChEMBL.issue"] != ""
    )
    start_match = (frame["PubMed.StartPage"] == frame["ChEMBL.first_page"]) & (
        frame["ChEMBL.first_page"] != ""
    )
    end_match = (frame["PubMed.EndPage"] == frame["ChEMBL.last_page"]) & (
        frame["ChEMBL.last_page"] != ""
    )

    journal_value = np.where(journal_match, frame["ChEMBL.journal"], "unknown")
    volume_value = np.where(volume_match, frame["ChEMBL.volume"], "unknown")
    issue_value = np.where(issue_match, frame["ChEMBL.issue"], "")
    start_value = np.where(start_match, frame["ChEMBL.first_page"], "unknown")
    end_value = np.where(end_match, frame["ChEMBL.last_page"], "")

    frame["reference"] = (
        pd.Series(journal_value, index=frame.index)
        + ", "
        + pd.Series(volume_value, index=frame.index)
        + "("
        + pd.Series(issue_value, index=frame.index)
        + "), p."
        + pd.Series(start_value, index=frame.index)
        + "-"
        + pd.Series(end_value, index=frame.index)
    )

    frame["invalid.reference"] = frame["reference"].str.contains("unknown", na=False)

    year_completed = frame["PubMed.YearCompleted"].map(to_text)
    year_revised = frame["PubMed.YearRevised"].map(to_text)
    chembl_year = frame["ChEMBL.year"].map(to_text)
    frame["year"] = np.where(
        ~year_completed.map(null_or_empty),
        year_completed.map(pad4),
        np.where(
            ~year_revised.map(null_or_empty),
            year_revised.map(pad4),
            chembl_year.map(pad4),
        ),
    )

    month_completed = frame["PubMed.MonthCompleted"].map(to_text)
    month_revised = frame["PubMed.MonthRevised"].map(to_text)
    frame["month"] = np.where(
        ~month_completed.map(null_or_empty),
        month_completed.map(pad2),
        np.where(~month_revised.map(null_or_empty), month_revised.map(pad2), "00"),
    )

    day_completed = frame["PubMed.DayCompleted"].map(to_text)
    day_revised = frame["PubMed.DayRevised"].map(to_text)
    frame["day"] = np.where(
        ~day_completed.map(null_or_empty),
        day_completed.map(pad2),
        np.where(~day_revised.map(null_or_empty), day_revised.map(pad2), "00"),
    )

    frame["completed"] = (
        pd.Series(frame["year"], index=frame.index)
        + "-"
        + pd.Series(frame["month"], index=frame.index)
        + "-"
        + pd.Series(frame["day"], index=frame.index)
    )

    frame["ChEMBL.pubmed_id"] = frame["ChEMBL.pubmed_id"].map(pad_pmid8)
    frame["sortorder.document"] = (
        frame["PubMed.ISSN"].map(to_text)
        + ":"
        + frame["completed"]
        + ":"
        + frame["ChEMBL.pubmed_id"].map(to_text)
    )

    frame = frame.rename(columns={"PubMed.PMID": "PMID", "ChEMBL.doi": "doi"})
    frame["PMID"] = frame["PMID"].map(to_text)
    frame["doi"] = frame["doi"].map(to_text)

    frame["invalid"] = (
        frame["invalid.doi"] | frame["invalid.PMID"] | frame["invalid.reference"]
    )

    drop_candidates = [
        column for column in STAGE_REMOVED_COLUMNS if column in frame.columns
    ]
    frame = frame.drop(columns=drop_candidates)

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
    frame = frame.rename(columns=rename_map)
    frame["doi"] = frame["doi"].map(to_text)

    secondary_rename = {
        "chemical_list": "PubMed.chemical.list",
        "mesh.qualifiers": "PubMed.mesh.qualifiers",
        "mesh.descriptors": "PubMed.mesh.descriptors",
    }
    frame = frame.rename(columns=secondary_rename)

    for column in LOWERCASE_LIST_COLUMNS:
        if column in frame.columns:
            frame[column] = frame[column].map(safe_lower)

    frame = frame.merge(reference, on="document_chembl_id", how="left")

    review_values: list[Any] = []
    for current, doctype in zip(frame["review"], frame["doctype_review"], strict=False):
        current_value = None if pd.isna(current) else bool(current)
        doctype_value = None if pd.isna(doctype) else bool(doctype)
        if current_value is True or doctype_value is True:
            review_values.append(True)
        elif current_value is False and doctype_value is False:
            review_values.append(False)
        else:
            review_values.append(pd.NA)
    frame["review"] = pd.Series(review_values, index=frame.index, dtype="boolean")

    experimental_values: list[Any] = []
    for value in frame["review"]:
        if value is pd.NA or pd.isna(value):
            experimental_values.append(pd.NA)
        else:
            experimental_values.append(not bool(value))
    frame["experimental"] = pd.Series(
        experimental_values, index=frame.index, dtype="boolean"
    )

    frame = frame.drop(columns=["is_experimental_doc"], errors="ignore")

    missing = [column for column in FINAL_COLUMN_ORDER if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing expected columns after harmonisation: {missing}")

    frame = frame.loc[:, FINAL_COLUMN_ORDER]
    if frame.columns.has_duplicates:
        frame = frame.loc[:, ~frame.columns.duplicated()]
    frame = frame.sort_values("completed", ascending=True, kind="mergesort")
    frame.reset_index(drop=True, inplace=True)

    return frame


# ===== Convenience wrapper ==================================================


def postprocess_file(
    input_path: Path | str,
    output_path: Path | str,
    *,
    cfg: IoCfg,
    ref_document_path: Path | str | None = None,
) -> Path:
    """Read a CSV, apply :func:`postprocess_documents` and write the result."""

    input_path = Path(input_path)
    output_path = Path(output_path)

    if output_path.is_dir() or output_path.suffix == "":
        destination = output_path / f"{OUTPUT_PREFIX}{input_path.name}"
    else:
        name = output_path.name
        if not name.startswith(OUTPUT_PREFIX):
            destination = output_path.with_name(f"{OUTPUT_PREFIX}{name}")
        else:
            destination = output_path

    frame = pd.read_csv(
        input_path,
        dtype=str,
        encoding=UTF8_ENCODING,
        sep=cfg.csv_sep,
        keep_default_na=True,
    )

    processed = postprocess_documents(
        frame,
        ref_document_path=ref_document_path,
    )

    write_csv_deterministic(
        processed,
        destination,
        col_order=list(processed.columns),
        key_cols=["completed", "document_chembl_id"],
        sep=cfg.csv_sep,
        encoding=UTF8_ENCODING,
        cfg=None,
    )

    return destination


__all__ = [
    "postprocess_documents",
    "postprocess_file",
    "normalize_journal",
    "null_or_empty",
    "pad4",
    "pad2",
    "pad_pmid8",
    "eq_text",
    "ne_text",
]
