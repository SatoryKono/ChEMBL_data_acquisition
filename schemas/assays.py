"""Schema definitions for assay data."""

from __future__ import annotations

from pandera import Check, Column, DataFrameSchema

AssaysSchema: DataFrameSchema = DataFrameSchema(
    {
        "assay_chembl_id": Column(str, required=True),
        "document_chembl_id": Column(str, required=True),
        "target_chembl_id": Column(str, required=False),
        "year": Column(int, Check.in_range(1900, 2100), required=True),
        "month": Column(int, Check.in_range(1, 12), required=True),
    }
)

"""pandera.DataFrameSchema: Validation schema for assays."""
