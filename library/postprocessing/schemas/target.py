from __future__ import annotations

import pandas as pd

from .common import TableSchema

TARGET_SCHEMA = TableSchema(
    name="TargetPostprocessSchema",
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
        "organism": "string",
        "target_class": "string",
        "protein_family": "string",
        "synonyms": "string",
        "pipeline_version": "string",
    },
    nullable_columns=(
        "organism",
        "target_class",
        "protein_family",
        "synonyms",
        "pipeline_version",
    ),
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
    sort_by=("target_chembl_id",),
)


def validate_targets(df: pd.DataFrame, *, context: str = "targets") -> pd.DataFrame:
    """Validate ``df`` against :data:`TARGET_SCHEMA`."""

    return TARGET_SCHEMA.validate(df, context=context)


__all__ = ["TARGET_SCHEMA", "validate_targets"]
