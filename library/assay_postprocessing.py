"""Assay table post-processing utilities.

This module provides helper functions to transform assay data retrieved from
ChEMBL.  The main feature adds an ``assay_with_same_target`` column that counts
how many assays share the same ``document_chembl_id`` and
``target_chembl_id`` combination.  The transformation mimics the behaviour of a
Power Query script used in the original workflow.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

import pandas as pd

logger = logging.getLogger(__name__)


def _validate_columns(df: pd.DataFrame, required: Iterable[str]) -> None:
    """Ensure that *df* contains all *required* columns.

    Parameters
    ----------
    df:
        DataFrame to validate.
    required:
        Iterable of mandatory column names.

    Raises
    ------
    ValueError
        If any of the required columns is missing from ``df``.
    """

    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"missing required columns: {', '.join(sorted(missing))}")


def postprocess_assays(df: pd.DataFrame) -> pd.DataFrame:
    """Augment assay data with per-target assay counts.

    The function adds a new column named ``assay_with_same_target`` to ``df``.  The
    value represents how many assays share the same ``document_chembl_id`` and
    ``target_chembl_id`` pair as the current row.

    Parameters
    ----------
    df:
        DataFrame produced by :func:`library.chembl_library.get_assays_all`.

    Returns
    -------
    pandas.DataFrame
        Copy of the input with the additional ``assay_with_same_target`` column.
    """

    _validate_columns(df, ["document_chembl_id", "target_chembl_id"])
    df = df.copy()
    counts = (
        df.groupby(["document_chembl_id", "target_chembl_id"])
        .size()
        .rename("assay_with_same_target")
    )
    logger.debug("Calculated counts for %d document/target groups", len(counts))
    df = df.merge(counts, on=["document_chembl_id", "target_chembl_id"], how="left")
    return df


def postprocess_file(
    input_path: Path | str,
    output_path: Path | str,
    *,
    sep: str = ",",
    encoding: str = "utf8",
) -> None:
    """Read an assay CSV, post-process and write the result.

    Parameters
    ----------
    input_path:
        Path to the CSV file produced by ``get_assay_data.py``.
    output_path:
        Destination path for the cleaned CSV file.
    sep:
        Field delimiter of the CSV files.
    encoding:
        Text encoding of the CSV files.
    """

    df = pd.read_csv(input_path, sep=sep, encoding=encoding, dtype=str)
    processed = postprocess_assays(df)
    processed.to_csv(output_path, index=False, sep=sep, encoding=encoding)
