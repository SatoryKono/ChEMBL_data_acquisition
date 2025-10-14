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


def merge_series_prefer_left(
    left: pd.Series[Any], right: pd.Series[Any]
) -> pd.Series[Any]:
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

    original_index = left.index

    left_work = left.copy()
    right_work = right.copy()

    requires_alignment = (
        not original_index.equals(right.index)
        or not left_work.index.is_unique
        or not right_work.index.is_unique
    )

    if requires_alignment:
        left_work = _ensure_unique_index(left_work, always_multi=True)
        right_work = _ensure_unique_index(right_work, always_multi=True)
        if not left_work.index.equals(right_work.index):
            right_work = right_work.reindex(left_work.index)

    result = left_work.copy()
    missing_mask = result.isna()
    if missing_mask.any():
        result.loc[missing_mask] = right_work.loc[missing_mask]

    if isinstance(result.index, pd.MultiIndex):
        result.index = original_index
    return result


def _ensure_unique_index(
    series: pd.Series[Any], *, always_multi: bool = False
) -> pd.Series[Any]:
    """Return ``series`` with a unique index preserving order.

    Pandas disallows :meth:`Series.reindex` when the axis contains duplicate
    labels. To align two series that may include duplicates we temporarily
    promote their index to a ``MultiIndex`` made of the original label and a
    per-label counter. The original order is preserved which allows a
    positionally stable alignment.
    """

    if series.index.is_unique and not always_multi:
        return series.copy()

    counts = series.groupby(level=0).cumcount()
    unique = series.copy()
    unique.index = pd.MultiIndex.from_arrays([series.index, counts])
    return unique


__all__ = ["json_normalize_pyarrow", "merge_series_prefer_left", "read_csv_pyarrow"]
