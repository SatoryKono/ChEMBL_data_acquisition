"""Writer utilities for CSV-based data pipelines."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from itertools import chain
from pathlib import Path

import pandas as pd

from ..common.csv_utils import write_csv_chunks_deterministic
from ..common.log import logger
from ..config import Config


def write_csv(
    df: pd.DataFrame | Iterable[pd.DataFrame],
    path: str | Path,
    *,
    cfg: Config,
    sep: str | None = None,
    encoding: str | None = None,
    key_cols: Iterable[str] | None = None,
    col_order: Iterable[str] | None = None,
    chunksize: int | None = None,
) -> Path:
    """Write ``df`` or dataframe chunks to ``path`` as CSV and store metadata."""

    sep = sep or cfg.io.csv_sep
    encoding = encoding or cfg.io.csv_encoding
    from ..common.csv_utils import write_csv_deterministic

    if isinstance(df, pd.DataFrame):
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

    iterator = iter(df)
    try:
        first = next(iterator)
    except StopIteration:
        empty = pd.DataFrame()
        return write_csv(
            empty,
            path,
            cfg=cfg,
            sep=sep,
            encoding=encoding,
            key_cols=key_cols,
            col_order=col_order,
            chunksize=chunksize,
        )

    if not isinstance(first, pd.DataFrame):
        raise TypeError("write_csv expects DataFrame or iterable of DataFrames")

    key_cols_list = list(key_cols) if key_cols is not None else sorted(first.columns)
    missing_keys = [col for col in key_cols_list if col not in first.columns]
    if missing_keys:
        logger.error(
            "missing_key_columns",
            requested=key_cols_list,
            missing=missing_keys,
        )
        raise ValueError(f"Missing key columns: {missing_keys}")

    if col_order is None:
        col_order_list = sorted(first.columns)
    else:
        col_order_list = list(col_order)

    chunk_iter: Iterator[pd.DataFrame] = chain([first], iterator)
    if chunksize is None:
        return write_csv_chunks_deterministic(
            chunk_iter,
            path,
            col_order=col_order_list,
            key_cols=key_cols_list,
            sep=sep,
            encoding=encoding,
            cfg=cfg,
        )

    return write_csv_chunks_deterministic(
        chunk_iter,
        path,
        col_order=col_order_list,
        key_cols=key_cols_list,
        chunksize=chunksize,
        sep=sep,
        encoding=encoding,
        cfg=cfg,
    )
