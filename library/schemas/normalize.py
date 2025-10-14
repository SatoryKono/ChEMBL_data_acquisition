"""Normalization helpers for schema dataframes.

This module provides utility functions for cleaning raw data before
validation against the :mod:`library.schemas` definitions.  The helpers ensure
consistent relation operators, unify units and trim identifier values.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import pandas as pd

_RELATION_MAP: dict[str, str] = {"<": "<=", ">": ">=", "=": "=", "<=": "<=", ">=": ">="}


def _apply_if_string(
    series: pd.Series[Any], func: Callable[[str], str]
) -> pd.Series[Any]:
    """Apply ``func`` to string elements of ``series``.

    Non-string values are returned unchanged.
    """
    if pd.api.types.is_bool_dtype(series):
        return series
    return series.map(lambda x: func(x) if isinstance(x, str) else x)


def _normalize_common(df: pd.DataFrame) -> pd.DataFrame:
    """Return a normalised copy of ``df``.

    The function standardises relation operators (``<`` → ``<=``, ``>`` → ``>=``),
    converts the micro sign ``μM`` to ``uM`` and strips whitespace from identifier
    columns.  Only string values are modified.
    """
    result = df.copy()
    for col in result.columns:
        series = result[col]
        series = _apply_if_string(series, lambda s: s.replace("μM", "uM"))
        if "relation" in col.lower():
            series = _apply_if_string(
                series, lambda s: _RELATION_MAP.get(s.strip(), s.strip())
            )
        if "id" in col.lower():
            # Ensure identifier columns are consistently treated as strings
            series = series.astype(pd.StringDtype())
            series = _apply_if_string(series, str.strip)
            series = series.astype(pd.StringDtype())
        result[col] = series
    return result


def normalize_activities(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise an activities dataframe.

    Parameters
    ----------
    df:
        Raw activities data.

    Returns
    -------
    pandas.DataFrame
        Normalised copy of ``df``.
    """
    return _normalize_common(df)


def normalize_assays(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise an assays dataframe.

    Parameters
    ----------
    df:
        Raw assays data.

    Returns
    -------
    pandas.DataFrame
        Normalised copy of ``df``.
    """
    return _normalize_common(df)


def normalize_documents(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise a documents dataframe.

    Parameters
    ----------
    df:
        Raw documents data.

    Returns
    -------
    pandas.DataFrame
        Normalised copy of ``df``.
    """
    return _normalize_common(df)


def normalize_cell_lines(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise a cell lines dataframe."""

    return _normalize_common(df)


def normalize_tissues(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise a tissues dataframe."""

    return _normalize_common(df)


def normalize_targets(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise a targets dataframe.

    Parameters
    ----------
    df:
        Raw targets data.

    Returns
    -------
    pandas.DataFrame
        Normalised copy of ``df``.
    """
    return _normalize_common(df)


def normalize_testitems(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise a test items dataframe.

    Parameters
    ----------
    df:
        Raw test items data.

    Returns
    -------
    pandas.DataFrame
        Normalised copy of ``df``.
    """
    return _normalize_common(df)


__all__ = [
    "normalize_activities",
    "normalize_assays",
    "normalize_documents",
    "normalize_cell_lines",
    "normalize_tissues",
    "normalize_targets",
    "normalize_testitems",
]
