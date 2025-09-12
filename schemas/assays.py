"""Schema definitions for assay data."""

from __future__ import annotations

import pandera.pandas as pa
from pandera import Check

AssaysSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "assay_chembl_id": pa.Column(str, required=True),
        "document_chembl_id": pa.Column(str, required=True),
        "target_chembl_id": pa.Column(str, required=True),
        "year": pa.Column(int, required=True, checks=Check.ge(0)),
        "month": pa.Column(int, required=True, checks=Check.between(1, 12)),
    }
)

"""pa.DataFrameSchema: Validation schema for assays."""
