"""Utilities for normalising multifunctional enzyme annotations."""

from __future__ import annotations

import pandas as pd


def normalise_multifunctional(series: pd.Series | None) -> pd.Series:
    """Return a nullable boolean series describing multifunctional enzymes."""

    if series is None:
        return pd.Series(pd.NA, index=pd.Index([], dtype="int64"), dtype="boolean")
    text = series.astype("string").str.strip().str.lower()
    mapped = text.map({"true": True, "false": False})
    return mapped.astype("boolean")
