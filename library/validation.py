"""Validation helpers for tabular datasets."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from .log import logger


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
    columns = list(required)
    logger.info("validate_start", extra={"stage": "validate_start", "columns": columns})
    missing = [col for col in columns if col not in df.columns]
    if missing:
        logger.info(
            "validate_done",
            extra={"stage": "validate_done", "columns": columns, "missing": missing},
        )
        raise ValueError(f"missing columns: {', '.join(missing)}")
    logger.info("validate_done", extra={"stage": "validate_done", "columns": columns})
