from __future__ import annotations

import pandas as pd

from .common import TableSchema

ACTIVITY_SCHEMA = TableSchema(
    name="ActivityPostprocessSchema",
    required_columns=(
        "activity_id",
        "molecule_chembl_id",
        "assay_chembl_id",
        "standard_type",
        "standard_relation",
        "standard_value",
        "standard_units",
    ),
    optional_columns=(
        "data_validity_comment",
        "activity_comment",
        "quality_flag",
        "pipeline_version",
    ),
    dtypes={
        "activity_id": "Int64",
        "molecule_chembl_id": "string",
        "assay_chembl_id": "string",
        "standard_type": "string",
        "standard_relation": "string",
        "standard_value": pd.Float64Dtype(),
        "standard_units": "string",
        "data_validity_comment": "string",
        "activity_comment": "string",
        "quality_flag": pd.BooleanDtype(),
        "pipeline_version": "string",
    },
    nullable_columns=(
        "data_validity_comment",
        "activity_comment",
        "quality_flag",
        "pipeline_version",
    ),
    column_order=(
        "activity_id",
        "molecule_chembl_id",
        "assay_chembl_id",
        "standard_type",
        "standard_relation",
        "standard_value",
        "standard_units",
        "data_validity_comment",
        "activity_comment",
        "quality_flag",
        "pipeline_version",
    ),
    sort_by=("molecule_chembl_id", "assay_chembl_id", "activity_id"),
)


def validate_activities(df: pd.DataFrame, *, context: str = "activities") -> pd.DataFrame:
    """Validate ``df`` against :data:`ACTIVITY_SCHEMA`."""

    return ACTIVITY_SCHEMA.validate(df, context=context)


__all__ = ["ACTIVITY_SCHEMA", "validate_activities"]
