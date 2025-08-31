"""Utilities for I/O operations.

This module contains helper functions to read input datasets and write
Parquet files in batches.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

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


def read_input(path: Path, sep: str, encoding: str) -> pd.DataFrame:
    """Read the input CSV/TSV file into a DataFrame.

    Parameters
    ----------
    path:
        Path to the input file.
    sep:
        Field separator.
    encoding:
        File encoding.

    Returns
    -------
    pd.DataFrame
        Loaded data with all columns as strings.
    """

    return pd.read_csv(path, sep=sep, encoding=encoding, dtype=str).fillna("")
