"""Schema specification for assay postprocessing."""

from __future__ import annotations

from library.postprocess.common import DataFrameSchema, validate_schema

ASSAY_SCHEMA = DataFrameSchema(
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
    },
    sort_by=("assay_chembl_id",),
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
)


def validate_assays(df, *, context: str = "assays"):
    """Validate ``df`` using :data:`ASSAY_SCHEMA`."""

    return validate_schema(df, ASSAY_SCHEMA, context=context)


__all__ = ["ASSAY_SCHEMA", "validate_assays"]
