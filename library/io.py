"""I/O helpers for CSV-based data pipelines.

This module centralises reading and writing of CSV files to ensure
consistent handling of delimiters, encodings and error reporting across
command line utilities.
"""

from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path
from typing import Iterable
import subprocess
import sys

import pandas as pd
import yaml

from . import validation
from .config import Config, IoCfg


def read_ids(
    path: str | Path,
    *,
    column: str,
    cfg: IoCfg,
    sep: str | None = None,
    encoding: str | None = None,
) -> list[str]:
    """Return identifier values from ``column`` in ``path``.

    Parameters
    ----------
    path:
        Location of the CSV file.
    column:
        Name of the column that contains identifiers.
    cfg:
        I/O configuration providing default CSV parameters.
    sep:
        Field delimiter used in the CSV file. Defaults to ``cfg.csv_sep``.
    encoding:
        Character encoding of the CSV file. Defaults to ``cfg.csv_encoding``.

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
    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding
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
    cfg: IoCfg,
    sep: str | None = None,
    encoding: str | None = None,
    required_columns: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Load a CSV file into a :class:`pandas.DataFrame` with optional schema validation.

    Parameters
    ----------
    path:
        Location of the CSV file.
    cfg:
        I/O configuration providing default CSV parameters.
    sep:
        Field delimiter used in the CSV file. Defaults to ``cfg.csv_sep``.
    encoding:
        Character encoding of the CSV file. Defaults to ``cfg.csv_encoding``.
    required_columns:
        Optional list of column names that must be present in the loaded
        DataFrame. A :class:`ValueError` is raised if any are missing.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the CSV contents.

    """
    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding
    df = pd.read_csv(path, sep=sep, encoding=encoding)
    if required_columns is not None:
        validation.validate_columns(df, required_columns)
    return df


def write_csv(
    df: pd.DataFrame,
    path: str | Path,
    *,
    cfg: Config,
    sep: str | None = None,
    encoding: str | None = None,
) -> None:
    """Write ``df`` to ``path`` as CSV and store metadata.

    Parameters
    ----------
    df:
        DataFrame to serialise.
    path:
        Destination file path.
    cfg:
        Full configuration used for metadata sidecars.
    sep:
        Field delimiter used in the CSV file. Defaults to ``cfg.io.csv_sep``.
    encoding:
        Character encoding of the CSV file. Defaults to ``cfg.io.csv_encoding``.

    """
    sep = sep or cfg.io.csv_sep
    encoding = encoding or cfg.io.csv_encoding
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, sep=sep, encoding=encoding)
    _write_meta(Path(path), cfg)


def default_output_path(input_path: str | Path) -> Path:
    """Return the default output path for ``input_path``.

    The generated name follows the pattern
    ``output_<stem>_YYYYMMDD.csv`` and is placed next to
    ``input_path``.
    """
    inp = Path(input_path)
    date_str = datetime.now().strftime("%Y%m%d")
    return inp.with_name(f"output_{inp.stem}_{date_str}.csv")


def _git_sha() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except Exception:  # pragma: no cover - git may be unavailable
        return "unknown"


def _write_meta(path: Path, cfg: Config) -> None:
    meta = {
        "git_sha": _git_sha(),
        "command": " ".join(sys.argv),
        "config": cfg.to_dict(),
    }
    meta_path = Path(str(path) + ".meta.yaml")
    with meta_path.open("w", encoding="utf8") as fh:
        yaml.safe_dump(meta, fh, sort_keys=False)
