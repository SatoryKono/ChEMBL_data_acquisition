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
from .config import Config, IoCfg, _serialize_paths
from .log import logger


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
    key_cols: Iterable[str] | None = None,
) -> None:
    """Write ``df`` to ``path`` as CSV and store metadata.

    Columns are sorted alphabetically and rows lexicographically to ensure
    deterministic output. Values are written using Unix line endings and a
    consistent floating-point representation.

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
    key_cols:
        Optional columns to determine row order. When provided, rows are
        sorted by these columns. Otherwise all columns are used, preserving
        the previous deterministic behaviour.

    """
    sep = sep or cfg.io.csv_sep
    encoding = encoding or cfg.io.csv_encoding
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    sorted_df = df.sort_index(axis=1)
    if key_cols is not None:
        # Only sort by the specified key columns for determinism.
        sorted_df = sorted_df.sort_values(by=list(key_cols))
    else:
        # Fall back to the original behaviour of sorting by all columns.
        sorted_df = sorted_df.sort_values(by=list(sorted_df.columns))
    sorted_df.to_csv(
        path,
        index=False,
        sep=sep,
        encoding=encoding,
        lineterminator="\n",
        float_format="%.6g",
    )
    _write_meta(Path(path), cfg)


def default_output_path(input_path: str | Path, cfg: IoCfg) -> Path:
    """Return the default output path for ``input_path``.

    Parameters
    ----------
    input_path:
        Source file path used only for deriving the stem name.
    cfg:
        I/O configuration providing the base directory for generated files.

    Returns
    -------
    pathlib.Path
        Path pointing inside ``cfg.output_dir`` following the
        ``output_<stem>_YYYYMMDD.csv`` naming scheme.
    """
    inp = Path(input_path)
    date_str = datetime.now().strftime("%Y%m%d")
    return Path(cfg.output_dir) / f"output_{inp.stem}_{date_str}.csv"


def _git_sha() -> str:
    """Return the current Git commit hash.

    The command is limited to a short timeout to avoid hanging when ``git``
    is unavailable. If the call times out or fails, ``"unknown"`` is
    returned and a warning is logged.
    """

    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
        return result.stdout.strip()
    except subprocess.TimeoutExpired as exc:
        logger.warning("git command timed out: %s", exc)
    except Exception as exc:  # pragma: no cover - git may be unavailable
        logger.warning("unable to determine git SHA: %s", exc)
    return "unknown"


def _write_meta(path: Path, cfg: Config) -> None:
    """Write YAML metadata alongside the output file.

    Paths inside the configuration are converted to plain strings so that the
    resulting YAML does not contain :class:`~pathlib.Path` objects.
    """

    meta = {
        "git_sha": _git_sha(),
        "command": " ".join(sys.argv),
        "config": _serialize_paths(cfg.to_dict()),
    }
    meta_path = Path(f"{path}.meta.yaml")
    with meta_path.open("w", encoding="utf8") as fh:
        yaml.safe_dump(meta, fh, sort_keys=False)
