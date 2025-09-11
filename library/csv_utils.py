"""CSV utilities for deterministic output.

This module provides :func:`write_csv_deterministic` which writes a
:pandas.DataFrame to disk in a reproducible manner.  Column order, row
sorting and value serialisation are normalised to guarantee that repeated
runs produce identical files.
"""

from __future__ import annotations

from datetime import date, datetime
import hashlib
import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence

import pandas as pd
import yaml
from pandas.api import types as ptypes

from .config import Config, _serialize_paths


logger = logging.getLogger(__name__)


def _git_sha() -> str:
    """Return the current Git commit hash or ``"unknown"`` if unavailable."""

    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
        return result.stdout.strip()
    except subprocess.TimeoutExpired:
        logger.warning("git rev-parse timed out")
        return "unknown"
    except Exception:  # pragma: no cover - git may be unavailable
        return "unknown"


def _write_meta(path: Path, cfg: Config | None) -> Path:
    """Write a sidecar YAML file with basic provenance information."""

    meta = {
        "git_sha": _git_sha(),
        "command": " ".join(sys.argv),
        "config": _serialize_paths(cfg.to_dict()) if cfg is not None else {},
    }
    meta_path = Path(f"{path}.meta.yaml")
    with meta_path.open("w", encoding="utf8") as fh:
        yaml.safe_dump(meta, fh, sort_keys=False)
    return meta_path


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
    sep: str = ",",
    encoding: str = "utf-8-sig",
    cfg: Config | None = None,
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

    Notes
    -----
    The input ``df`` is mutated in-place (column order, sorting and normalisation)
    to minimise memory usage. Callers that need to preserve the original data
    should pass ``df.copy()``.

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

    # Operations below mutate ``df`` directly to avoid unnecessary copies.
    # Callers should pass ``df.copy()`` if the original frame must remain
    # unchanged.
    work = df

    # Determine column order using ``reindex`` which can avoid materialising
    # a full copy when ``copy=False`` is passed.  This keeps memory overhead
    # proportional to the number of columns rather than the entire DataFrame.
    if col_order is not None:
        head = [c for c in col_order if c in work.columns]
        tail = sorted(c for c in work.columns if c not in head)
        work = work.reindex(columns=head + tail, copy=False)
    else:
        work = work.reindex(columns=sorted(work.columns), copy=False)

    # Sort rows deterministically in-place using a stable algorithm.
    sort_cols = list(key_cols) if key_cols is not None else list(work.columns)
    work.sort_values(by=sort_cols, kind="mergesort", inplace=True)

    # Normalise bool and date columns without creating intermediary frames.
    for col in work.columns:
        s = work[col]
        if ptypes.is_bool_dtype(s):
            work[col] = _normalise_bool(s)
        else:
            work[col] = _normalise_dates(s)

    # Write via temporary file for atomicity
    with tempfile.NamedTemporaryFile(
        "w", encoding=encoding, newline="\n", delete=False, dir=str(out_path.parent)
    ) as fh:
        tmp_path = Path(fh.name)
        work.to_csv(
            fh,
            index=False,
            float_format="%.6g",
            na_rep="",
            chunksize=chunksize,
            sep=sep,
        )
    os.replace(tmp_path, out_path)
    _write_meta(out_path, cfg)
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
