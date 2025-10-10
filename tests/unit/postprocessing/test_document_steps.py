from __future__ import annotations

import pandas as pd
import pytest

from library.postprocessing.documents.steps import (
    enrich_document_publication_year,
    finalize_document_records,
    normalize_document_fields,
)


@pytest.mark.unit
def test_normalize_document_fields__applies_whitespace_and_unicode() -> None:
    frame = pd.DataFrame(
        {
            "Title": ["  Cafe\u0301  "],
            "Journal": ["  Journal   of   Testing  "],
            "Doc_Type": ["  PREPRINT  "],
        }
    )

    result = normalize_document_fields(
        frame,
        trim_whitespace=True,
        normalise_unicode=True,
    )

    assert list(result.columns) == ["title", "journal", "doc_type"]
    assert result.loc[0, "title"] == "Café"
    assert result.loc[0, "journal"] == "Journal of Testing"
    assert result.loc[0, "doc_type"] == "PREPRINT"
    assert frame.loc[0, "Title"] == "  Cafe\u0301  "


@pytest.mark.unit
def test_normalize_document_fields__skips_whitespace_when_disabled() -> None:
    frame = pd.DataFrame({"Title": ["  Example  "]})

    result = normalize_document_fields(frame, trim_whitespace=False)

    assert result.loc[0, "title"] == "  Example  "


@pytest.mark.unit
def test_enrich_document_publication_year__uses_year_and_fallback() -> None:
    frame = pd.DataFrame({"year": ["2010", pd.NA, ""]})

    result = enrich_document_publication_year(frame, fallback_year=1999)

    expected = pd.Series([2010, 1999, 1999], dtype="Int64")
    pd.testing.assert_series_equal(
        result["publication_year"], expected, check_names=False
    )


@pytest.mark.unit
def test_enrich_document_publication_year__prefers_external_sources_when_requested() -> (
    None
):
    frame = pd.DataFrame(
        {
            "year": ["2001", "2002"],
            "pubmed.yearcompleted": ["2005", pd.NA],
        }
    )

    result = enrich_document_publication_year(frame, prefer_doi_year=True)

    expected = pd.Series([2005, 2002], dtype="Int64")
    pd.testing.assert_series_equal(
        result["publication_year"], expected, check_names=False
    )


@pytest.mark.unit
def test_enrich_document_publication_year__validates_fallback_range() -> None:
    frame = pd.DataFrame({"year": [pd.NA]})

    with pytest.raises(ValueError, match="supported range"):
        enrich_document_publication_year(frame, fallback_year=1400)


@pytest.mark.unit
def test_finalize_document_records__adds_missing_required_columns() -> None:
    frame = pd.DataFrame({"document_chembl_id": ["CHEMBL1"]})

    result = finalize_document_records(frame)

    assert list(result.columns[:3]) == [
        "document_chembl_id",
        "title",
        "doc_type",
    ]
    assert result.loc[0, "document_chembl_id"] == "CHEMBL1"
    assert pd.isna(result.loc[0, "title"])
    assert pd.isna(result.loc[0, "doc_type"])
    assert result["document_chembl_id"].dtype == "string"
    assert result["title"].dtype == "string"
    assert result["doc_type"].dtype == "string"


@pytest.mark.unit
def test_finalize_document_records__deduplicates_when_requested() -> None:
    frame = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL1", "CHEMBL1", "CHEMBL2"],
            "title": ["A", "B", "C"],
            "doc_type": ["Journal", "Preprint", "Journal"],
        }
    )

    result = finalize_document_records(frame, ensure_unique_ids=True)

    assert result["document_chembl_id"].tolist() == ["CHEMBL1", "CHEMBL2"]
    assert result["title"].tolist() == ["A", "C"]


@pytest.mark.unit
def test_finalize_document_records__fills_missing_identifier_from_prefixed_column() -> (
    None
):
    frame = pd.DataFrame(
        {
            "chembl.document_chembl_id": ["CHEMBL42", "CHEMBL99"],
            "title": ["First", "Second"],
        }
    )

    result = finalize_document_records(frame, ensure_unique_ids=True)

    assert result["document_chembl_id"].tolist() == ["CHEMBL42", "CHEMBL99"]
    assert result["chembl.document_chembl_id"].tolist() == ["CHEMBL42", "CHEMBL99"]


@pytest.mark.unit
def test_finalize_document_records__skips_validation_when_disabled() -> None:
    frame = pd.DataFrame({"unexpected": [1]})

    result = finalize_document_records(frame, enforce_schema=False)

    assert "unexpected" in result.columns
    # Missing required columns are injected even when validation is disabled.
    assert all(
        column in result.columns
        for column in ("document_chembl_id", "title", "doc_type")
    )
