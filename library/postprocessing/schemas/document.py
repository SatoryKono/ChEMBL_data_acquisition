from __future__ import annotations

import pandas as pd

from .common import TableSchema

DOCUMENT_SCHEMA = TableSchema(
    name="DocumentPostprocessSchema",
    required_columns=(
        "document_chembl_id",
        "title",
        "doc_type",
    ),
    optional_columns=(
        "year",
        "doi",
        "journal",
        "abstract",
        "publication_year",
        "pipeline_version",
    ),
    dtypes={
        "document_chembl_id": "string",
        "title": "string",
        "doc_type": "string",
        "year": "string",
        "doi": "string",
        "journal": "string",
        "abstract": "string",
        "publication_year": "Int64",
        "pipeline_version": "string",
    },
    nullable_columns=(
        "year",
        "doi",
        "journal",
        "abstract",
        "publication_year",
        "pipeline_version",
    ),
    column_order=(
        "document_chembl_id",
        "title",
        "doc_type",
        "year",
        "publication_year",
        "doi",
        "journal",
        "abstract",
        "pipeline_version",
    ),
    sort_by=("document_chembl_id",),
)


def validate_documents(df: pd.DataFrame, *, context: str = "documents") -> pd.DataFrame:
    """Validate ``df`` against :data:`DOCUMENT_SCHEMA`."""

    return DOCUMENT_SCHEMA.validate(df, context=context)


__all__ = ["DOCUMENT_SCHEMA", "validate_documents"]
