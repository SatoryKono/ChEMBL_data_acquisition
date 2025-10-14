"""Pandera schema describing activity postprocessing outputs."""

from __future__ import annotations

from library._compat.pandera import pa

from .common import boolean_column, float_column, int_column, string_column

__all__ = ["ActivitySchema"]


ActivitySchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "activity_id": int_column(required=True, nullable=False),
        "molecule_chembl_id": string_column(required=True, nullable=False),
        "assay_chembl_id": string_column(required=True, nullable=False),
        "standard_type": string_column(required=True, nullable=False),
        "standard_relation": string_column(required=True, nullable=False),
        "standard_value": float_column(required=True, nullable=True),
        "standard_units": string_column(required=True, nullable=False),
        "data_validity_comment": string_column(required=False, nullable=True),
        "activity_comment": string_column(required=False, nullable=True),
        "quality_flag": boolean_column(required=False, nullable=True),
        "pipeline_version": string_column(required=False, nullable=True),
    }
)
