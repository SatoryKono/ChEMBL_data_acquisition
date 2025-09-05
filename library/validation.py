"""Validation helpers for tabular datasets."""

from __future__ import annotations

from typing import Iterable

import pandas as pd


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
