"""Schema specification for document postprocessing."""

from __future__ import annotations

from library.postprocess.common import DataFrameSchema, validate_schema

DOCUMENT_SCHEMA = DataFrameSchema(
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
        "publication_year": "Int64",
    },
    sort_by=("document_chembl_id",),
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
)


def validate_documents(df, *, context: str = "documents"):
    """Validate ``df`` using :data:`DOCUMENT_SCHEMA`."""

    return validate_schema(df, DOCUMENT_SCHEMA, context=context)


__all__ = ["DOCUMENT_SCHEMA", "validate_documents"]
