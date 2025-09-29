"""Writer utilities for CSV-based data pipelines."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pandas as pd

from ..config import Config
from ..log import logger


def write_csv(
    df: pd.DataFrame,
    path: str | Path,
    *,
    cfg: Config,
    sep: str | None = None,
    encoding: str | None = None,
    key_cols: Iterable[str] | None = None,
    col_order: Iterable[str] | None = None,
    chunksize: int | None = None,
) -> Path:
    """Write ``df`` to ``path`` as CSV and store metadata."""

    sep = sep or cfg.io.csv_sep
    encoding = encoding or cfg.io.csv_encoding
    key_cols_list = list(key_cols) if key_cols is not None else sorted(df.columns)
    missing_keys = [col for col in key_cols_list if col not in df.columns]
    if missing_keys:
        logger.error(
            "missing_key_columns",
            requested=key_cols_list,
            missing=missing_keys,
        )
        raise ValueError(f"Missing key columns: {missing_keys}")
    col_order_list = list(col_order) if col_order is not None else None
    from ..csv_utils import write_csv_deterministic

    return write_csv_deterministic(
        df,
        path,
        key_cols=key_cols_list,
        col_order=col_order_list,
        chunksize=chunksize,
        sep=sep,
        encoding=encoding,
        cfg=cfg,
    )
