"""Schema specification for target postprocessing."""
from __future__ import annotations

from library.postprocessing.common import DataFrameSchema, validate_schema

TARGET_SCHEMA = DataFrameSchema(
    required_columns=(
        "target_chembl_id",
        "pref_name",
        "target_type",
    ),
    optional_columns=(
        "organism",
        "target_class",
        "protein_family",
        "synonyms",
        "pipeline_version",
    ),
    dtypes={
        "target_chembl_id": "string",
        "pref_name": "string",
        "target_type": "string",
    },
    sort_by=("target_chembl_id",),
    column_order=(
        "target_chembl_id",
        "pref_name",
        "target_type",
        "organism",
        "target_class",
        "protein_family",
        "synonyms",
        "pipeline_version",
    ),
)


def validate_targets(df, *, context: str = "targets"):
    """Validate ``df`` using :data:`TARGET_SCHEMA`."""

    return validate_schema(df, TARGET_SCHEMA, context=context)


__all__ = ["TARGET_SCHEMA", "validate_targets"]
