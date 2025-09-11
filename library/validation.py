"""Validation helpers for tabular datasets."""

from __future__ import annotations

from typing import Iterable, Mapping

import pandas as pd
import pandas.api.types as ptypes

from .log import logger


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


def validate_schema(df: pd.DataFrame, schema: Mapping[str, str]) -> None:
    """Validate presence and dtype of DataFrame columns.

    Parameters
    ----------
    df:
        DataFrame to inspect.
    schema:
        Mapping of column names to expected pandas dtypes (as strings).

    Raises
    ------
    ValueError
        If required columns are missing.
    TypeError
        If a column has an unexpected dtype.

    """
    schema_dict = dict(schema)
    logger.info(
        "validate_start", extra={"stage": "validate_start", "schema": schema_dict}
    )
    validate_columns(df, schema_dict.keys())
    for column, dtype in schema_dict.items():
        if not ptypes.is_dtype_equal(df[column].dtype, dtype):
            logger.info(
                "validate_done",
                extra={
                    "stage": "validate_done",
                    "schema": schema_dict,
                    "column": column,
                    "expected": dtype,
                    "actual": str(df[column].dtype),
                },
            )
            raise TypeError(
                f"column {column!r} has dtype {df[column].dtype}, expected {dtype}"
            )
    logger.info(
        "validate_done", extra={"stage": "validate_done", "schema": schema_dict}
    )
