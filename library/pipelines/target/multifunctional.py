"""Helpers for inferring the ``multifunctional_enzyme`` flag."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import pandas as pd

DEFAULT_TRUE_TOKENS: frozenset[str] = frozenset(
    {
        "multifunctional",
        "multi-functional",
        "multi functional",
    }
)
DEFAULT_FALSE_TOKENS: frozenset[str] = frozenset({"", "-", "nan", "none"})
DEFAULT_COLUMNS: tuple[str, ...] = (
    "protein_classifications",
    "IUPHAR_class",
    "IUPHAR_subclass",
)

__all__ = [
    "DEFAULT_COLUMNS",
    "DEFAULT_FALSE_TOKENS",
    "DEFAULT_TRUE_TOKENS",
    "append_multifunctional_flag",
    "infer_multifunctional",
]


def _normalise(value: Any) -> str:
    if isinstance(value, pd.Series):  # defensive
        raise TypeError("_normalise expects a scalar value")
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    text = str(value).strip()
    lowered = text.lower()
    return lowered


def _token_matches(token: str, true_tokens: Iterable[str]) -> bool:
    for candidate in true_tokens:
        if candidate and candidate in token:
            return True
    return False


def infer_multifunctional(
    record: Mapping[str, Any],
    *,
    columns: Sequence[str] = DEFAULT_COLUMNS,
    true_tokens: Iterable[str] = DEFAULT_TRUE_TOKENS,
    false_tokens: Iterable[str] = DEFAULT_FALSE_TOKENS,
) -> bool:
    """Infer the multifunctional flag from a mapping of annotations."""

    for column in columns:
        value = _normalise(record.get(column))
        if not value or value in false_tokens:
            continue
        if _token_matches(value, true_tokens):
            return True
    return False


def append_multifunctional_flag(
    df: pd.DataFrame,
    *,
    output_column: str = "multifunctional_enzyme",
    columns: Sequence[str] = DEFAULT_COLUMNS,
    true_tokens: Iterable[str] = DEFAULT_TRUE_TOKENS,
    false_tokens: Iterable[str] = DEFAULT_FALSE_TOKENS,
) -> pd.DataFrame:
    """Return ``df`` with a derived ``multifunctional_enzyme`` column."""

    result = df.copy()
    if result.empty:
        result[output_column] = pd.Series(dtype="boolean")
        return result

    def _infer(row: pd.Series) -> bool:
        return infer_multifunctional(
            row,
            columns=columns,
            true_tokens=true_tokens,
            false_tokens=false_tokens,
        )

    inferred = result.apply(_infer, axis=1)
    result[output_column] = inferred.astype("boolean")
    return result
