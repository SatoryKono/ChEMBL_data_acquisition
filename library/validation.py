"""Validation helpers for tabular datasets."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Protocol

import pandas as pd
from pandera.errors import SchemaErrors

from .log import logger

__all__ = [
    "ValidationResult",
    "validate_columns",
    "validate_activities",
    "validate_assays",
    "validate_testitems",
]

_FAILURE_COLUMNS = [
    "schema_context",
    "column",
    "check",
    "check_number",
    "failure_case",
    "index",
]


@dataclass(slots=True)
class ValidationResult:
    """Outcome of applying a pandera schema to a dataframe."""

    data: pd.DataFrame
    failure_cases: pd.DataFrame
    schema_name: str

    @property
    def has_failures(self) -> bool:
        """Return ``True`` when validation detected failing rows."""

        return not self.failure_cases.empty


def _empty_failure_cases() -> pd.DataFrame:
    """Return an empty dataframe matching pandera failure case layout."""

    return pd.DataFrame(columns=_FAILURE_COLUMNS)


def _unique_indices(values: pd.Series) -> list[Any]:
    """Return cleaned list of failing indices from ``values``."""

    if values.empty:
        return []
    cleaned: list[Any] = []
    for value in pd.unique(values.dropna()):
        if isinstance(value, float) and value.is_integer():
            cleaned.append(int(value))
        else:
            cleaned.append(value)
    return cleaned


def _combine_failure_cases(cases: list[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate collected failure case frames preserving columns."""

    if not cases:
        return _empty_failure_cases()
    return pd.concat(cases, ignore_index=True, sort=False).reindex(
        columns=_FAILURE_COLUMNS, fill_value=pd.NA
    )


class _SupportsDataFrameValidate(Protocol):
    """Protocol describing the subset of pandera schema API we rely on."""

    def validate(
        self, dataframe: pd.DataFrame, *args: Any, **kwargs: Any
    ) -> pd.DataFrame:
        """Validate *dataframe* and return the cleansed result."""


def _validate_with_schema(
    df: pd.DataFrame,
    *,
    schema: _SupportsDataFrameValidate,
    schema_name: str,
    return_result: bool,
) -> ValidationResult | pd.DataFrame:
    """Apply ``schema`` to ``df`` dropping failing rows when possible."""

    logger.info(
        "schema_validate_start",
        extra={"schema": schema_name, "rows": len(df)},
    )
    current = df
    collected: list[pd.DataFrame] = []
    while True:
        try:
            validated = schema.validate(current, lazy=True)
        except SchemaErrors as exc:
            failure_cases = exc.failure_cases.copy()
            indices = _unique_indices(
                failure_cases["index"]
                if "index" in failure_cases
                else pd.Series(dtype=object)
            )
            if not indices:
                logger.info(
                    "schema_validate_error",
                    extra={
                        "schema": schema_name,
                        "rows": len(current),
                        "failures": len(failure_cases),
                    },
                )
                raise
            collected.append(failure_cases)
            current = current.drop(index=indices, errors="ignore")
            continue
        collected_cases = _combine_failure_cases(collected)
        dropped = len(df) - len(validated)
        logger.info(
            "schema_validate_done",
            extra={
                "schema": schema_name,
                "rows": len(validated),
                "dropped": dropped,
                "failures": len(collected_cases),
            },
        )
        result = ValidationResult(validated, collected_cases, schema_name)
        return result if return_result else result.data


def validate_columns(df: pd.DataFrame, required: Iterable[str]) -> None:
    """Ensure that all ``required`` columns exist in ``df``.

    Parameters
    ----------
    df:
        DataFrame to inspect.
    required:
        Names of columns that must be present in ``df``.

    Raises
    ------
    ValueError
        If any columns are missing.

    """
    columns = list(required)
    logger.info("validate_start", extra={"stage": "validate_start", "columns": columns})
    missing = [col for col in columns if col not in df.columns]
    if missing:
        logger.info(
            "validate_done",
            extra={"stage": "validate_done", "columns": columns, "missing": missing},
        )
        raise ValueError(f"missing columns: {', '.join(missing)}")
    logger.info("validate_done", extra={"stage": "validate_done", "columns": columns})


def validate_activities(
    df: pd.DataFrame, *, return_result: bool = False
) -> ValidationResult | pd.DataFrame:
    """Validate activities dataframe using :data:`schemas.ActivitiesSchema`."""

    from schemas import ActivitiesSchema

    return _validate_with_schema(
        df,
        schema=ActivitiesSchema,
        schema_name="ActivitiesSchema",
        return_result=return_result,
    )


def validate_assays(
    df: pd.DataFrame, *, return_result: bool = False
) -> ValidationResult | pd.DataFrame:
    """Validate assay dataframe using :data:`schemas.AssaysSchema`."""

    from schemas import AssaysSchema

    return _validate_with_schema(
        df,
        schema=AssaysSchema,
        schema_name="AssaysSchema",
        return_result=return_result,
    )


def validate_testitems(
    df: pd.DataFrame, *, return_result: bool = False
) -> ValidationResult | pd.DataFrame:
    """Validate test item dataframe using :data:`schemas.TestitemsSchema`."""

    from schemas import TestitemsSchema

    return _validate_with_schema(
        df,
        schema=TestitemsSchema,
        schema_name="TestitemsSchema",
        return_result=return_result,
    )
