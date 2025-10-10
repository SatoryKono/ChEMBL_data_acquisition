"""Shared text conversion helpers used across pipelines and QA checks."""

from __future__ import annotations

from typing import Any, cast

import numpy as np
import pandas as pd

UTF8_ENCODING = "utf-8"


def to_text(value: Any, *, encoding: str = UTF8_ENCODING) -> str:
    """Convert *value* to a string while mapping null-like values to ``""``."""

    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, bytes):
        return value.decode(encoding, errors="ignore")
    if isinstance(value, np.bool_ | bool):
        return "true" if bool(value) else "false"
    if isinstance(value, np.integer | int):
        return str(int(value))
    if isinstance(value, np.floating | float):
        numeric = float(value)
        if numeric != numeric:  # NaN check
            return ""
        if numeric.is_integer():
            return str(int(numeric))
        return str(numeric)
    if pd.isna(cast(object, value)):
        return ""
    return str(value)


__all__ = ["to_text", "UTF8_ENCODING"]
