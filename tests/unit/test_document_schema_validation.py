"""Unit tests for the document Pandera schema helpers."""

from __future__ import annotations

import pandas as pd
import pytest
from pandera.errors import SchemaError

from library.schemas.document_schema import DocumentSchema, validate_document_frame


@pytest.mark.unit
def test_validate_document_frame__accepts_base_columns_only() -> None:
    frame = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL1"],
            "doi": ["10.1000/example"],
            "pubmed_id": [12345],
            "title": ["Example"],
            "doc_type": ["journal"],
            "journal": ["Journal"],
            "year": [2024],
            "doi_key": ["10.1000/example"],
        }
    )

    validated = validate_document_frame(frame)

    assert validated.loc[0, "pubmed_id"] == "12345"
    assert str(validated.dtypes["pubmed_id"]) == "string"
    assert set(DocumentSchema.columns) >= set(validated.columns)


@pytest.mark.unit
def test_validate_document_frame__accepts_enrichment_columns() -> None:
    frame = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL2"],
            "doi": ["10.2000/test"],
            "pubmed_id": [pd.NA],
            "title": [pd.NA],
            "doc_type": [pd.NA],
            "journal": [pd.NA],
            "year": [pd.NA],
            "doi_key": ["10.2000/test"],
            "crossref_title": ["CrossRef"],
            "crossref_doc_type": ["journal-article"],
            "crossref_subject": ["Biology"],
            "crossref_error": [pd.NA],
            "openalex_title": ["OpenAlex"],
            "openalex_doc_type": ["article"],
            "openalex_type_crossref": ["journal-article"],
            "openalex_publication_year": [2023],
            "openalex_error": [pd.NA],
        }
    )

    validated = validate_document_frame(frame)

    assert validated.loc[0, "openalex_publication_year"] == "2023"
    assert str(validated.dtypes["crossref_title"]) == "string"


@pytest.mark.unit
def test_validate_document_frame__missing_required_column() -> None:
    frame = pd.DataFrame(
        {
            "doi": ["10.3000/example"],
            "pubmed_id": ["123"],
            "title": ["Example"],
            "doc_type": ["journal"],
            "journal": ["Journal"],
            "year": ["2024"],
            "doi_key": ["10.3000/example"],
        }
    )

    with pytest.raises(SchemaError):
        validate_document_frame(frame)
