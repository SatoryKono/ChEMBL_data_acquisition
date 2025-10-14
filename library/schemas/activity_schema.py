"""Pandera schema enforcing the normalized activity dataset contract."""

from __future__ import annotations

import pandas as pd

from library._compat.pandera import pa

from .common import float_column, int_column, string_column

__all__ = ["ACTIVITY_SCHEMA", "validate_activity_records", "ACTIVITY_COLUMN_ORDER"]


ACTIVITY_COLUMN_ORDER: tuple[str, ...] = (
    "activity_id",
    "assay_chembl_id",
    "target_chembl_id",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
)


ACTIVITY_SCHEMA: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "activity_id": int_column(required=True, nullable=False),
        "assay_chembl_id": string_column(required=True, nullable=False),
        "target_chembl_id": string_column(required=True, nullable=True),
        "standard_type": string_column(required=True, nullable=True),
        "standard_relation": string_column(required=True, nullable=True),
        "standard_value": float_column(required=True, nullable=True),
        "standard_units": string_column(required=True, nullable=True),
    },
    strict=True,
    coerce=True,
)


def validate_activity_records(df: pd.DataFrame) -> pd.DataFrame:
    """Validate *df* against :data:`ACTIVITY_SCHEMA` returning the coerced frame."""

    return ACTIVITY_SCHEMA.validate(df, lazy=True)
