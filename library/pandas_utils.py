"""Pandas compatibility helpers.

These utility functions provide fallbacks for features that are only available in
newer versions of pandas (>=2.0), such as the ``dtype_backend`` argument.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, cast

import pandas as pd


def json_normalize_pyarrow(
    data: Sequence[Mapping[str, Any]] | Mapping[str, Any],
    *,
    dtype_backend: str = "pyarrow",
    **kwargs: Any,
) -> pd.DataFrame:
    """Normalize semi-structured JSON data.

    This wrapper calls :func:`pandas.json_normalize` with the ``dtype_backend``
    parameter when supported. If running on an older pandas version that does
    not accept ``dtype_backend``, the argument is omitted to maintain backward
    compatibility.

    Parameters
    ----------
    data:
        Nested data to be normalized.
    dtype_backend:
        Backend to use for dtypes, passed to :func:`pandas.json_normalize` when
        available. Defaults to ``"pyarrow"``.
    **kwargs:
        Additional keyword arguments forwarded to
        :func:`pandas.json_normalize`.

    Returns
    -------
    pandas.DataFrame
        Normalized data.
    """

    json_normalize = cast(Any, pd.json_normalize)
    try:
        # pandas >= 2.0 supports ``dtype_backend``
        return cast(
            pd.DataFrame,
            json_normalize(
                data,
                dtype_backend=dtype_backend,
                **kwargs,
            ),
        )
    except TypeError:
        # For pandas < 2.0, fall back to the default behaviour
        return cast(pd.DataFrame, json_normalize(data, **kwargs))


def read_csv_pyarrow(
    *args: Any,
    dtype_backend: str = "pyarrow",
    **kwargs: Any,
) -> pd.DataFrame:
    """Read a CSV file with optional ``dtype_backend`` support.

    Parameters
    ----------
    *args:
        Positional arguments forwarded to :func:`pandas.read_csv`.
    dtype_backend:
        Backend to use for dtypes, passed to :func:`pandas.read_csv` when
        available. Defaults to ``"pyarrow"``.
    **kwargs:
        Additional keyword arguments forwarded to :func:`pandas.read_csv`.

    Returns
    -------
    pandas.DataFrame
        Parsed CSV data.
    """

    read_csv = cast(Any, pd.read_csv)
    try:
        return cast(
            pd.DataFrame,
            read_csv(
                *args,
                dtype_backend=dtype_backend,
                **kwargs,
            ),
        )
    except (TypeError, ImportError):
        return cast(pd.DataFrame, read_csv(*args, **kwargs))


def merge_series_prefer_left(left: pd.Series, right: pd.Series) -> pd.Series:
    """Return ``left`` with missing entries populated from ``right``.

    The function mirrors :meth:`pandas.Series.combine_first` for the subset of
    behaviour relied upon in the codebase while avoiding the deprecated
    concatenation of empty arrays that triggers a ``FutureWarning`` in pandas
    2.2+. Only missing values (``NaN``/``NA``) are replaced, preserving the
    existing handling of empty strings for object columns.

    Parameters
    ----------
    left:
        Primary series whose non-missing values take precedence.
    right:
        Fallback series used to fill missing entries in ``left``.

    Returns
    -------
    pandas.Series
        A copy of ``left`` with missing positions filled from ``right`` where
        possible.
    """

    if left.empty and right.empty:
        return left.copy()

    if not left.index.equals(right.index):
        right = right.reindex(left.index)

    result = left.copy()
    missing_mask = result.isna()
    if missing_mask.any():
        result.loc[missing_mask] = right.loc[missing_mask]
    return result


__all__ = ["json_normalize_pyarrow", "merge_series_prefer_left", "read_csv_pyarrow"]
