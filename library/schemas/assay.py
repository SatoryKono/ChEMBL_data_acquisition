"""Pandera schema describing assay postprocessing outputs."""

from __future__ import annotations

from library._compat.pandera import pa

from .common import boolean_column, float_column, string_column

__all__ = ["AssaySchema"]


AssaySchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "assay_chembl_id": string_column(required=True, nullable=False),
        "assay_type": string_column(required=True, nullable=False),
        "assay_test_type": string_column(required=True, nullable=False),
        "description": string_column(required=True, nullable=True),
        "assay_format": string_column(required=False, nullable=True),
        "target_chembl_id": string_column(required=False, nullable=True),
        "confidence_score": float_column(required=False, nullable=True),
        "is_confirmatory": boolean_column(required=False, nullable=True),
        "pipeline_version": string_column(required=False, nullable=True),
    }
)
