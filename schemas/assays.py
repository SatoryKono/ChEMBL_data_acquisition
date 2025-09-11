"""Schema definitions for assay data."""

from __future__ import annotations

import pandera.pandas as pa

AssaysSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "assay_chembl_id": pa.Column(str, required=True),
        "document_chembl_id": pa.Column(str, required=True),
        "target_chembl_id": pa.Column(str, required=False),
        "year": pa.Column(int, pa.Check.in_range(1900, 2100), required=True),
        "month": pa.Column(int, pa.Check.in_range(1, 12), required=True),
    }
)

"""pa.DataFrameSchema: Validation schema for assays."""
