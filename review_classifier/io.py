from __future__ import annotations

"""I/O utilities for document classification."""

from pathlib import Path
from typing import Iterator

import pandas as pd


def read_documents(path: Path, chunksize: int | None = None) -> pd.DataFrame | Iterator[pd.DataFrame]:
    """Read the main documents table.

    Parameters
    ----------
    path:
        Path to ``all_documents.csv``.
    chunksize:
        Optional chunk size for streaming.

    Returns
    -------
    pandas.DataFrame or Iterator[pandas.DataFrame]
        Loaded dataframe or iterator over chunks.
    """
    if not Path(path).exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    return pd.read_csv(path, chunksize=chunksize)


def read_mesh_probabilities(path: Path) -> dict[str, float]:
    """Load experimental MeSH probabilities into a mapping.

    Parameters
    ----------
    path:
        CSV with columns ``Custom`` and ``experimental_probability``.

    Returns
    -------
    dict[str, float]
        Mapping ``term -> probability``.
    """
    if not Path(path).exists():
        raise FileNotFoundError(f"Mesh probability file not found: {path}")
    df = pd.read_csv(path)
    required = {"Custom", "experimental_probability"}
    if not required.issubset(df.columns):
        raise ValueError("experimental_mesh.csv must contain columns 'Custom' and 'experimental_probability'")
    df["Custom"] = df["Custom"].astype(str).str.lower().str.strip()
    return dict(zip(df["Custom"], df["experimental_probability"].astype(float)))


def write_output(df: pd.DataFrame, path: Path, fmt: str = "csv") -> None:
    """Write the classified dataframe to disk.

    Parameters
    ----------
    df:
        DataFrame to save.
    path:
        Output path.
    fmt:
        ``csv`` or ``parquet``.
    """
    if fmt == "csv":
        df.to_csv(path, index=False)
    elif fmt == "parquet":
        df.to_parquet(path, index=False)
    else:
        raise ValueError(f"Unsupported output format: {fmt}")
