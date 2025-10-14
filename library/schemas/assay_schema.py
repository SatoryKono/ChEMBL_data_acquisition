"""Pandera schema describing assay exports produced by the pipeline."""

from __future__ import annotations

import pandas as pd

from library._compat.pandera import pa

from .common import int_column, string_column

__all__ = ["AssayDataSchema"]


AssayDataSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "assay_chembl_id": string_column(required=True, nullable=False),
        "assay_type": string_column(required=False, nullable=True),
        "assay_test_type": string_column(required=False, nullable=True),
        "target_chembl_id": string_column(required=False, nullable=True),
        "created_on": string_column(required=False, nullable=True),
        "updated_on": string_column(required=False, nullable=True),
        "assay_strain": string_column(required=False, nullable=True),
        "assay_group": string_column(required=False, nullable=True),
        "accession": string_column(required=False, nullable=True),
        "timestamp_utc": pa.Column(
            pd.DatetimeTZDtype(tz="UTC"),
            required=False,
            nullable=True,
            coerce=True,
        ),
        "year": int_column(required=False, nullable=True),
    }
)
