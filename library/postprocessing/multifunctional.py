"""Utilities for normalising multifunctional enzyme annotations."""

from __future__ import annotations

from typing import Any

import pandas as pd


def normalise_multifunctional(
    series: pd.Series[Any] | None,
) -> pd.Series[bool]:
    """Return a nullable boolean series describing multifunctional enzymes."""

    if series is None:
        index = pd.Index([], dtype="int64")
        data = [pd.NA] * len(index)
        return pd.Series(data, index=index, dtype="boolean")
    text = series.astype("string").str.strip().str.lower()
    mapped = text.map({"true": True, "false": False})
    return mapped.astype("boolean")
