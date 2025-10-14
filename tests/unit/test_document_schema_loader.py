"""Tests for :mod:`library.schemas.document_schema`."""

from __future__ import annotations

import pandas as pd
import pytest

from library.schemas.document_schema import DOCUMENT_SCHEMA, DocumentSchema, validate_document_frame


@pytest.mark.unit
def test_document_schema_instance_is_pandera_schema() -> None:
    assert isinstance(DOCUMENT_SCHEMA, DocumentSchema)
    assert DOCUMENT_SCHEMA.columns


@pytest.mark.unit
def test_validate_document_frame_preserves_columns() -> None:
    frame = pd.DataFrame(
        [
            {
                "document_chembl_id": "DOC1",
                "title": "title",
                "abstract": "",
                "doi": "10.1000/example",
                "year": "2024",
                "journal": "Nature",
                "journal_abbrev": "Nat",
                "volume": "10",
                "issue": "1",
                "first_page": "1",
                "last_page": "12",
                "pubmed_id": "1234",
                "extra_column": "value",
            }
        ]
    )

    validated = validate_document_frame(frame)

    assert "extra_column" in validated.columns
    assert validated.loc[0, "document_chembl_id"] == "DOC1"
