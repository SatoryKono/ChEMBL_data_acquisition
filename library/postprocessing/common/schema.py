from __future__ import annotations
from typing import Any, cast

import pandas as pd

from library.postprocessing.schemas.common import TableSchema, validate_with_pandera

from .types import SchemaValidationError

DataFrameSchema = TableSchema


def validate_schema(
    df: pd.DataFrame, schema: DataFrameSchema, *, context: str
) -> pd.DataFrame:
    """Validate ``df`` using the Pandera-backed ``schema``."""

    validated = validate_with_pandera(df, schema, context=context)

    ordered = validated.copy(deep=True)

    if schema.column_order:
        ordered_columns = [
            column for column in schema.column_order if column in ordered.columns
        ]
        remaining = sorted(
            column for column in ordered.columns if column not in ordered_columns
        )
        ordered = ordered.loc[:, [*ordered_columns, *remaining]].copy(deep=True)

    if schema.sort_by:
        sort_columns = [column for column in schema.sort_by if column in ordered.columns]
        if sort_columns:
            ordered = ordered.sort_values(
                sort_columns,
                kind="mergesort",
            ).reset_index(drop=True)

    return ordered


def coerce_types(df: pd.DataFrame, schema: DataFrameSchema) -> pd.DataFrame:
    """Attempt to coerce DataFrame columns according to ``schema.dtypes``."""

    coerced = df.copy(deep=True)
    for column, expected_type in schema.dtypes.items():
        if column not in coerced.columns:
            continue
        try:
            if isinstance(expected_type, str):
                coerced[column] = coerced[column].astype(cast(str, expected_type))
            else:
                coerced[column] = coerced[column].astype(cast(Any, expected_type))
        except (TypeError, ValueError) as exc:
            raise SchemaValidationError(
                column,
                f"failed to coerce to {expected_type!r}: {exc}",
                cause=exc,
            ) from exc
    return coerced


__all__ = ["DataFrameSchema", "coerce_types", "validate_schema"]
