from __future__ import annotations

import pandas as pd

from .common import TableSchema

ASSAY_SCHEMA = TableSchema(
    name="AssayPostprocessSchema",
    required_columns=(
        "assay_chembl_id",
        "assay_type",
        "assay_test_type",
        "description",
    ),
    optional_columns=(
        "assay_format",
        "target_chembl_id",
        "confidence_score",
        "is_confirmatory",
        "pipeline_version",
    ),
    dtypes={
        "assay_chembl_id": "string",
        "assay_type": "string",
        "assay_test_type": "string",
        "description": "string",
        "assay_format": "string",
        "target_chembl_id": "string",
        "confidence_score": pd.Float64Dtype(),
        "is_confirmatory": pd.BooleanDtype(),
        "pipeline_version": "string",
    },
    column_order=(
        "assay_chembl_id",
        "assay_type",
        "assay_test_type",
        "assay_format",
        "description",
        "target_chembl_id",
        "confidence_score",
        "is_confirmatory",
        "pipeline_version",
    ),
    sort_by=("assay_chembl_id",),
)


def validate_assays(df: pd.DataFrame, *, context: str = "assays") -> pd.DataFrame:
    """Validate ``df`` against :data:`ASSAY_SCHEMA`."""

    return ASSAY_SCHEMA.validate(df, context=context)


__all__ = ["ASSAY_SCHEMA", "validate_assays"]
