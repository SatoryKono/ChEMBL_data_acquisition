"""Schema definitions for activity postprocessing outputs."""

from __future__ import annotations

from library.postprocessing.pipeline.common import DataFrameSchema, validate_schema

ACTIVITY_SCHEMA = DataFrameSchema(
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
        "standard_units": "string",
    },
    sort_by=("molecule_chembl_id", "assay_chembl_id", "activity_id"),
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
)


def validate_activities(df, *, context: str = "activities"):
    """Validate ``df`` using :data:`ACTIVITY_SCHEMA`."""

    return validate_schema(df, ACTIVITY_SCHEMA, context=context)


__all__ = ["ACTIVITY_SCHEMA", "validate_activities"]
