"""CSV utilities for deterministic output.

This module provides :func:`write_csv_deterministic` which writes a
:pandas.DataFrame to disk in a reproducible manner.  Column order, row
sorting and value serialisation are normalised to guarantee that repeated
runs produce identical files.
"""

from __future__ import annotations

from datetime import date, datetime
import hashlib
import os
import tempfile
from pathlib import Path
from typing import Sequence

import pandas as pd
from pandas.api import types as ptypes


def _normalise_bool(series: pd.Series) -> pd.Series:
    """Return ``series`` with booleans converted to ``"true"``/``"false"``.

    ``pandas`` has two boolean dtypes: the classic ``bool`` which cannot
    represent missing values and the nullable ``boolean`` extension type. We
    map both explicit values to lowercase strings and leave missing values as
    ``<NA>`` which :func:`pandas.DataFrame.to_csv` will serialise using
    ``na_rep``.
    """

    return series.map({True: "true", False: "false"}).astype("string")


def _normalise_dates(series: pd.Series) -> pd.Series:
    """Return ``series`` with dates formatted as ``YYYY-MM-DD``.

    Both native ``datetime`` columns and object columns containing
    :class:`datetime.date` or :class:`datetime.datetime` values are
    supported. Non-date values are left untouched.
    """

    if ptypes.is_datetime64_any_dtype(series):
        return series.dt.strftime("%Y-%m-%d")

    if (
        ptypes.is_object_dtype(series)
        and series.dropna().map(lambda x: isinstance(x, (date, datetime))).all()
    ):
        return pd.to_datetime(series).dt.strftime("%Y-%m-%d")

    return series


def write_csv_deterministic(
    df: pd.DataFrame,
    path: str | Path,
    *,
    col_order: Sequence[str] | None = None,
    key_cols: Sequence[str] | None = None,
    chunksize: int | None = None,
) -> Path:
    """Serialise ``df`` to ``path`` as a deterministic CSV file.

    Parameters
    ----------
    df:
        DataFrame to be written.
    path:
        Destination file path.
    col_order:
        Optional preferred column ordering. Columns not listed here are
        appended in lexicographical order.
    key_cols:
        Optional sequence of column names defining row ordering.  If omitted
        all columns are used as sort keys.
    chunksize:
        Optional number of rows to write per chunk. Passing a value enables
        streaming output via :meth:`pandas.DataFrame.to_csv`, reducing peak
        memory usage for large tables. ``None`` (the default) writes the
        dataframe in a single call.

    Returns
    -------
    pathlib.Path
        Path to the written CSV file.
    """

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Validate optional chunksize to avoid infinite loops or crashes
    if chunksize is not None and chunksize <= 0:
        msg = f"chunksize must be a positive integer, got {chunksize}"
        raise ValueError(msg)

    work = df.copy()

    # Determine column order
    if col_order is not None:
        head = [c for c in col_order if c in work.columns]
        tail = sorted(c for c in work.columns if c not in head)
        work = work[head + tail]
    else:
        work = work[sorted(work.columns)]

    # Sort rows deterministically
    sort_cols = list(key_cols) if key_cols is not None else list(work.columns)
    work = work.sort_values(by=sort_cols, kind="mergesort")

    # Normalise bool and date columns
    for col in work.columns:
        s = work[col]
        if ptypes.is_bool_dtype(s):
            work[col] = _normalise_bool(s)
        else:
            work[col] = _normalise_dates(s)

    # Write via temporary file for atomicity
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8-sig", newline="\n", delete=False, dir=str(out_path.parent)
    ) as fh:
        tmp_path = Path(fh.name)
        work.to_csv(
            fh,
            index=False,
            float_format="%.6g",
            na_rep="",
            chunksize=chunksize,
        )
    os.replace(tmp_path, out_path)
    return out_path


def sha256_file(path: Path, *, block_size: int = 64 * 1024) -> str:
    """Return SHA-256 hash of ``path``.

    The file is read in chunks to avoid loading large files entirely into
    memory.

    Parameters
    ----------
    path:
        Path to the input file.
    block_size:
        Number of bytes read per iteration. Defaults to ``64 * 1024``.

    Returns
    -------
    str
        Hexadecimal SHA-256 digest of the file contents.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    """

    if not path.is_file():
        raise FileNotFoundError(f"File not found: {path}")

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(block_size), b""):
            digest.update(chunk)
    return digest.hexdigest()
