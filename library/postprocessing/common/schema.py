"""Utilities for declaratively defining and validating DataFrame schemas."""

from __future__ import annotations

import logging
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, cast

import pandas as pd
from pandera.errors import SchemaErrors

from library._compat.pandera import pa

from .types import SchemaValidationError

_LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class DataFrameSchema:
    """Declarative schema for postprocessing outputs."""

    required_columns: Sequence[str] = field(default_factory=tuple)
    optional_columns: Sequence[str] = field(default_factory=tuple)
    dtypes: Mapping[str, str | type] = field(default_factory=dict)
    sort_by: Sequence[str] | None = None
    column_order: Sequence[str] | None = None
    pandera_schema: pa.DataFrameSchema | None = None

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


def _summarise_failure_cases(cases: pd.DataFrame) -> str:
    """Return a short textual summary of pandera failure cases."""

    if cases.empty:
        return "no failure cases recorded"
    preview = cases.head(5).to_dict("records")
    return f"{len(cases)} failure(s): {preview}"


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

    current = df.copy(deep=True)

    if schema.pandera_schema is not None:
        schema_name = (
            getattr(schema.pandera_schema, "name", None)
            or schema.pandera_schema.__class__.__name__
        )
        _LOGGER.info(
            "Schema validation started",
            extra={"context": context, "schema": schema_name, "rows": len(current)},
        )
        try:
            current = schema.pandera_schema.validate(current, lazy=True)
        except SchemaErrors as exc:
            failure_cases = exc.failure_cases.copy()
            _LOGGER.error(
                "Schema validation failed",
                extra={
                    "context": context,
                    "schema": schema_name,
                    "error_count": len(failure_cases),
                    "errors_preview": failure_cases.head(5).to_dict("records"),
                },
            )
            message = _summarise_failure_cases(failure_cases)
            raise SchemaValidationError(context, message, cause=exc) from exc
        else:
            _LOGGER.info(
                "Schema validation succeeded",
                extra={
                    "context": context,
                    "schema": schema_name,
                    "rows": len(current),
                },
            )

    if schema.dtypes:
        mismatched: list[str] = []
        for column, expected_type in schema.dtypes.items():
            if column not in current.columns:
                continue
            dtype = current[column].dtype
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

    ordered_df = current.copy(deep=True)

    if schema.column_order:
        ordered_columns = [
            col for col in schema.column_order if col in ordered_df.columns
        ]
        remaining = [
            col for col in ordered_df.columns if col not in ordered_columns
        ]
        if ordered_columns:
            ordered_df = ordered_df.loc[:, ordered_columns + remaining]
        elif remaining:
            ordered_df = ordered_df.loc[:, remaining]

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
                dtype_name = cast(str, expected_type)
                coerced[column] = coerced[column].astype(dtype_name)
            else:
                dtype_object = cast(Any, expected_type)
                coerced[column] = coerced[column].astype(dtype_object)
        except (TypeError, ValueError) as exc:
            raise SchemaValidationError(
                column,
                f"failed to coerce to {expected_type!r}: {exc}",
                cause=exc,
            ) from exc
    return coerced


__all__ = ["DataFrameSchema", "coerce_types", "validate_schema"]
