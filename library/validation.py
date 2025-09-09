"""Validation helpers for tabular datasets."""

from __future__ import annotations

from typing import Iterable, Mapping

import pandas as pd
import pandas.api.types as ptypes


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
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"missing columns: {', '.join(missing)}")


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

    validate_columns(df, schema.keys())
    for column, dtype in schema.items():
        if not ptypes.is_dtype_equal(df[column].dtype, dtype):
            raise TypeError(
                f"column {column!r} has dtype {df[column].dtype}, expected {dtype}"
            )
