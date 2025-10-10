"""Utilities for declaratively defining and validating DataFrame schemas."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import pandas as pd

from .types import SchemaValidationError


@dataclass(frozen=True)
class DataFrameSchema:
    """Declarative schema for postprocessing outputs."""

    required_columns: Sequence[str] = field(default_factory=tuple)
    optional_columns: Sequence[str] = field(default_factory=tuple)
    dtypes: Mapping[str, str | type] = field(default_factory=dict)
    sort_by: Sequence[str] | None = None
    column_order: Sequence[str] | None = None

    def all_columns(self) -> Sequence[str]:
        ordered = list(self.required_columns) + list(self.optional_columns)
        if self.column_order:
            # Column order might include columns not listed above (derived fields).
            seen: set[str] = set()
            ordered = []
            for col in self.column_order:
                if col not in seen:
                    ordered.append(col)
                    seen.add(col)
        return ordered


def _missing_columns(df: pd.DataFrame, expected: Iterable[str]) -> list[str]:
    return [col for col in expected if col not in df.columns]


def validate_schema(
    df: pd.DataFrame, schema: DataFrameSchema, *, context: str
) -> pd.DataFrame:
    """Validate ``df`` against ``schema`` and return a defensive copy."""

    missing_required = _missing_columns(df, schema.required_columns)
    if missing_required:
        raise SchemaValidationError(
            context,
            f"missing required columns: {', '.join(sorted(missing_required))}",
        )

    if schema.dtypes:
        mismatched: list[str] = []
        for column, expected_type in schema.dtypes.items():
            if column not in df.columns:
                continue
            dtype = df[column].dtype
            if isinstance(expected_type, str):
                expected_name = expected_type
                if dtype.name != expected_type:
                    mismatched.append(
                        f"{column} (expected {expected_name}, got {dtype})"
                    )
            else:
                expected_dtype = pd.api.types.pandas_dtype(expected_type)
                expected_name = str(expected_dtype)
                if not pd.api.types.is_dtype_equal(dtype, expected_dtype):
                    mismatched.append(
                        f"{column} (expected {expected_name}, got {dtype})"
                    )
        if mismatched:
            raise SchemaValidationError(
                context,
                "; ".join(mismatched),
            )

    ordered_df = df.copy(deep=True)

    if schema.column_order:
        ordered_columns = [
            col for col in schema.column_order if col in ordered_df.columns
        ]
        remaining = sorted(
            col for col in ordered_df.columns if col not in ordered_columns
        )
        ordered_df = ordered_df.loc[:, ordered_columns + remaining]

    if schema.sort_by:
        sort_cols = [col for col in schema.sort_by if col in ordered_df.columns]
        if sort_cols:
            ordered_df = ordered_df.sort_values(
                sort_cols, kind="mergesort"
            ).reset_index(drop=True)

    return ordered_df


def coerce_types(df: pd.DataFrame, schema: DataFrameSchema) -> pd.DataFrame:
    """Attempt to coerce columns according to ``schema.dtypes``."""

    coerced = df.copy(deep=True)
    for column, expected_type in schema.dtypes.items():
        if column not in coerced.columns:
            continue
        try:
            if isinstance(expected_type, str):
                coerced[column] = coerced[column].astype(expected_type)
            else:
                coerced[column] = coerced[column].astype(expected_type)
        except (TypeError, ValueError) as exc:
            raise SchemaValidationError(
                column,
                f"failed to coerce to {expected_type!r}: {exc}",
                cause=exc,
            ) from exc
    return coerced


__all__ = ["DataFrameSchema", "coerce_types", "validate_schema"]
