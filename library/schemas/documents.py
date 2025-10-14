"""Schema definitions for document data."""

from __future__ import annotations

from library._compat.pandera import DataFrameSchema

from .document_spec import DOCUMENT_DECLARATION, DOCUMENT_SCHEMA_COLUMNS

DocumentsSchema: DataFrameSchema = DOCUMENT_DECLARATION.schema
"""DataFrameSchema: Validation schema for documents."""

__all__ = ["DocumentsSchema", "DOCUMENT_SCHEMA_COLUMNS"]
