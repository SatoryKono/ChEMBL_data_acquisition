"""Schema definitions for assay data."""

from __future__ import annotations

import pandera.pandas as pa

AssaysSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "assay_chembl_id": pa.Column(str, required=True),
        "document_chembl_id": pa.Column(str, required=True),
        "target_chembl_id": pa.Column(str, required=False),
        "year": pa.Column(int, required=False),
        "month": pa.Column(int, pa.Check.in_range(1, 12), required=False),
    }
)

"""pa.DataFrameSchema: Validation schema for assays."""
