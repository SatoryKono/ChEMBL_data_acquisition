"""Schema definitions for document data."""

from __future__ import annotations

from pandera import Check, Column, DataFrameSchema

DocumentsSchema: DataFrameSchema = DataFrameSchema(  # type: ignore[no-untyped-call]
    {
        "document_chembl_id": Column(str, required=True),
        "doi": Column(str, required=False),
        "title": Column(str, required=True),
        "year": Column(int, Check.in_range(1900, 2100), required=True),
        "month": Column(int, Check.in_range(1, 12), required=True),
        "day": Column(int, Check.in_range(1, 31), required=False),
        "citation": Column(int, Check.ge(0), required=False),
    }
)

"""pandera.DataFrameSchema: Validation schema for documents."""
