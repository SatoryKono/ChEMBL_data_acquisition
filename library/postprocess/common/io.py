"""I/O helpers for converting objects into :class:`pandas.DataFrame` instances."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import pandas as pd


def ensure_dataframe(data: object, *, copy: bool = False) -> pd.DataFrame:
    """Ensure that ``data`` is a :class:`~pandas.DataFrame` instance.

    Parameters
    ----------
    data:
        Any object that can be converted to a DataFrame.
    copy:
        Whether to return a defensive copy of the input DataFrame.
    """

    if isinstance(data, pd.DataFrame):
        return data.copy(deep=True) if copy else data

    if isinstance(data, Mapping):
        return pd.DataFrame(data)

    if isinstance(data, (Sequence, Iterable)):
        return pd.DataFrame(list(data))

    raise TypeError(f"Unsupported input for ensure_dataframe: {type(data)!r}")


def clone_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Return a deep copy of ``df`` ensuring metadata is preserved."""

    clone = df.copy(deep=True)
    clone.attrs.update(df.attrs)
    return clone


__all__ = ["clone_dataframe", "ensure_dataframe"]
