"""Schema definitions for document data."""

from __future__ import annotations

from .document_schema import DOCUMENT_SCHEMA, DOCUMENT_SCHEMA_COLUMNS as _ORDERED_COLUMNS
from .document_schema import DocumentSchema

DocumentsSchema: DocumentSchema = DOCUMENT_SCHEMA
"""DataFrame schema used to validate document metadata tables."""

DOCUMENT_SCHEMA_COLUMNS = _ORDERED_COLUMNS
"""Ordered column names provided by :class:`~library.schemas.document_schema.DocumentSchema`."""

__all__ = ["DocumentsSchema", "DOCUMENT_SCHEMA_COLUMNS"]
