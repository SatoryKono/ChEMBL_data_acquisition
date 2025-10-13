"""Schema definitions for document data."""

from __future__ import annotations

from library._compat.pandera import pa

from .document_spec import DOCUMENT_DECLARATION, DOCUMENT_SCHEMA_COLUMNS

DocumentsSchema: pa.DataFrameSchema = DOCUMENT_DECLARATION.schema
"""pa.DataFrameSchema: Validation schema for documents."""

__all__ = ["DocumentsSchema", "DOCUMENT_SCHEMA_COLUMNS"]
