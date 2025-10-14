"""Pandera schema for the normalized activity export."""

from __future__ import annotations

import pandas as pd

from library._compat.pandera import pa

from .common import float_column, int_column, string_column

__all__ = ["ACTIVITY_SCHEMA", "validate_activity"]


ACTIVITY_SCHEMA = pa.DataFrameSchema(
    {
        "activity_id": int_column(required=True, nullable=False),
        "assay_chembl_id": string_column(required=True, nullable=True),
        "target_chembl_id": string_column(required=True, nullable=True),
        "standard_type": string_column(required=True, nullable=True),
        "standard_relation": string_column(required=True, nullable=True),
        "standard_value": float_column(required=True, nullable=True),
        "standard_units": string_column(required=True, nullable=True),
    },
    coerce=True,
    strict=True,
)


def validate_activity(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate ``frame`` against :data:`ACTIVITY_SCHEMA`."""

    return ACTIVITY_SCHEMA.validate(frame, lazy=True)

