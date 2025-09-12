"""Pandas compatibility helpers.

These utility functions provide fallbacks for features that are only available in
newer versions of pandas (>=2.0), such as the ``dtype_backend`` argument.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import pandas as pd  # type: ignore[import-untyped]


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

    try:
        # pandas >= 2.0 supports ``dtype_backend``
        return pd.json_normalize(data, dtype_backend=dtype_backend, **kwargs)
    except TypeError:
        # For pandas < 2.0, fall back to the default behaviour
        return pd.json_normalize(data, **kwargs)


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

    try:
        return pd.read_csv(*args, dtype_backend=dtype_backend, **kwargs)
    except (TypeError, ImportError):
        return pd.read_csv(*args, **kwargs)


__all__ = ["json_normalize_pyarrow", "read_csv_pyarrow"]
