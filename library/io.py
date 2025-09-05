"""I/O helpers for CSV-based data pipelines.

This module centralises reading and writing of CSV files to ensure
consistent handling of delimiters, encodings and error reporting across
command line utilities.
"""

from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path

import pandas as pd


def read_ids(
    path: str | Path,
    *,
    column: str,
    sep: str = ",",
    encoding: str = "utf8",
) -> list[str]:
    """Return identifier values from ``column`` in ``path``.

    Parameters
    ----------
    path:
        Location of the CSV file.
    column:
        Name of the column that contains identifiers.
    sep:
        Field delimiter used in the CSV file. Defaults to ``","``.
    encoding:
        Character encoding of the CSV file. Defaults to ``"utf8"``.

    Returns
    -------
    list[str]
        Identifier values in the order they appear. Empty strings and
        ``"#N/A"`` markers are discarded.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the CSV file is malformed or ``column`` is missing.
    """
    try:
        with Path(path).open("r", encoding=encoding, newline="") as fh:
            reader = csv.DictReader(fh, delimiter=sep)
            if reader.fieldnames is None or column not in reader.fieldnames:
                raise ValueError(f"column '{column}' not found in {path}")
            ids: list[str] = []
            for row in reader:
                value = (row.get(column) or "").strip()
                if value and value != "#N/A":
                    ids.append(value)
            return ids
    except FileNotFoundError:
        raise
    except csv.Error as exc:
        raise ValueError(f"malformed CSV in file: {path}: {exc}") from exc


def read_csv(
    path: str | Path,
    *,
    sep: str = ",",
    encoding: str = "utf8",
) -> pd.DataFrame:
    """Load a CSV file into a :class:`pandas.DataFrame`.

    Parameters
    ----------
    path:
        Location of the CSV file.
    sep:
        Field delimiter used in the CSV file. Defaults to ``","``.
    encoding:
        Character encoding of the CSV file. Defaults to ``"utf8"``.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the CSV contents.
    """
    return pd.read_csv(path, sep=sep, encoding=encoding)


def write_csv(
    df: pd.DataFrame,
    path: str | Path,
    *,
    sep: str = ",",
    encoding: str = "utf8",
) -> None:
    """Write ``df`` to ``path`` as CSV.

    Parameters
    ----------
    df:
        DataFrame to serialise.
    path:
        Destination file path.
    sep:
        Field delimiter used in the CSV file. Defaults to ``","``.
    encoding:
        Character encoding of the CSV file. Defaults to ``"utf8"``.
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, sep=sep, encoding=encoding)


def default_output_path(input_path: str | Path) -> Path:
    """Return the default output path for ``input_path``.

    The generated name follows the pattern
    ``output_<stem>_YYYYMMDD.csv`` and is placed next to
    ``input_path``.
    """
    inp = Path(input_path)
    date_str = datetime.now().strftime("%Y%m%d")
    return inp.with_name(f"output_{inp.stem}_{date_str}.csv")
