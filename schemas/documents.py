"""Schema definitions for document data."""

from __future__ import annotations

import pandera.pandas as pa

DocumentsSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "document_chembl_id": pa.Column(str, required=True),
        "doi": pa.Column(str, required=False),
        "title": pa.Column(str, required=True),
        "year": pa.Column(int, pa.Check.in_range(1900, 2100), required=True),
        "month": pa.Column(int, pa.Check.in_range(1, 12), required=True),
        "day": pa.Column(int, pa.Check.in_range(1, 31), required=False),
        "citation": pa.Column(int, pa.Check.ge(0), required=False),
    }
)

"""pa.DataFrameSchema: Validation schema for documents."""
