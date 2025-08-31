"""Utilities for I/O operations.

This module contains helper functions to read input datasets and write
Parquet files in batches.
"""
from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Optional

import chardet
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


@dataclass
class ParquetBatchWriter:
    """Write pandas DataFrames to Parquet in batches.

    Parameters
    ----------
    path: Path
        Destination Parquet file.
    schema: Optional[pa.Schema]
        Schema to use for the Parquet file. If ``None`` the schema of the
        first written batch is used.
    """

    path: Path
    schema: Optional[pa.Schema] = None
    _writer: Optional[pq.ParquetWriter] = None

    def write(self, df: pd.DataFrame) -> None:
        """Append a batch to the Parquet file.

        Parameters
        ----------
        df:
            Data to append. The schema of the DataFrame must match the
            schema used for previous batches.
        """

        table = pa.Table.from_pandas(df, schema=self.schema, preserve_index=False)
        if self._writer is None:
            self.schema = table.schema
            self._writer = pq.ParquetWriter(self.path, self.schema)
        self._writer.write_table(table)

    def close(self) -> None:
        """Close the underlying writer."""

        if self._writer is not None:
            self._writer.close()
            self._writer = None


def _detect_encoding(path: Path) -> str:
    """Detect the encoding of ``path`` using :mod:`chardet`.

    Parameters
    ----------
    path:
        File to analyse.

    Returns
    -------
    str
        Detected encoding or ``"utf-8"`` as a fallback.
    """

    with path.open("rb") as fh:
        raw = fh.read(100_000)
    result = chardet.detect(raw)
    encoding = result.get("encoding") or "utf-8"
    logging.debug(
        "Detected encoding %s for %s (confidence %.2f)",
        encoding,
        path,
        result.get("confidence", 0.0),
    )
    return encoding


def read_input(path: Path, sep: str, encoding: Optional[str] = None) -> pd.DataFrame:
    """Read the input CSV/TSV file into a :class:`~pandas.DataFrame`.

    Parameters
    ----------
    path:
        Path to the input file.
    sep:
        Field separator.
    encoding:
        File encoding. If ``None`` the encoding is auto-detected.

    Returns
    -------
    pd.DataFrame
        Loaded data with all columns as strings.

    Raises
    ------
    ValueError
        If the file cannot be decoded with the detected/specified encoding.
    """

    if encoding is None:
        encoding = _detect_encoding(path)
    try:
        return pd.read_csv(path, sep=sep, encoding=encoding, dtype=str).fillna("")
    except UnicodeDecodeError as exc:
        raise ValueError(
            f"Could not decode {path} with encoding {encoding}. "
            "Pass a valid encoding via --encoding."
        ) from exc
