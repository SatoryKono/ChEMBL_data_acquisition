"""Tests validating CSV rows against DataFrame schemas."""

from __future__ import annotations

import pandas as pd
import pytest
from pandera.errors import SchemaError

from schemas.activities import ActivitiesSchema
from schemas.documents import DocumentsSchema


def test_activities_schema_enforces_float_dtype(
    activities_sample_df: pd.DataFrame,
) -> None:
    """`ActivitiesSchema` enforces column data types."""

    ActivitiesSchema.validate(activities_sample_df)

    bad = activities_sample_df.astype({"standard_value": str})
    with pytest.raises(SchemaError):
        ActivitiesSchema.validate(bad)


def test_documents_schema_requires_identifier(
    documents_sample_df: pd.DataFrame,
) -> None:
    """`DocumentsSchema` requires the document identifier column."""

    DocumentsSchema.validate(documents_sample_df)

    missing = documents_sample_df.drop(columns=["document_chembl_id"])
    with pytest.raises(SchemaError):
        DocumentsSchema.validate(missing)
