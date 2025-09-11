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

import pandas as pd

from .config import IoCfg
from .validation import validate_schema

logger = logging.getLogger(__name__)


def postprocess_assays(df: pd.DataFrame) -> pd.DataFrame:
    """Augment assay data with per-target assay counts.

    The function adds a new column named ``assay_with_same_target`` to ``df``.  The
    value represents how many assays share the same ``document_chembl_id`` and
    ``target_chembl_id`` pair as the current row.

    Parameters
    ----------
    df:
        DataFrame produced by :func:`library.chembl_library.get_assays`.

    Returns
    -------
    pandas.DataFrame
        Copy of the input with the additional ``assay_with_same_target`` column.

    """
    validate_schema(
        df,
        {
            "document_chembl_id": "object",
            "target_chembl_id": "object",
        },
    )
    group_cols = ["document_chembl_id", "target_chembl_id"]
    groups = df.groupby(group_cols)
    logger.debug("Calculated counts for %d document/target groups", groups.ngroups)
    result = df.copy()
    result["assay_with_same_target"] = groups["document_chembl_id"].transform("size")
    return result


def postprocess_file(
    input_path: Path | str,
    output_path: Path | str,
    *,
    cfg: IoCfg,
    sep: str | None = None,
    encoding: str | None = None,
) -> None:
    """Read an assay CSV, post-process and write the result.

    Parameters
    ----------
    input_path:
        Path to the CSV file produced by ``get_assay_data.py``.
    output_path:
        Destination path for the cleaned CSV file.
    cfg:
        I/O configuration providing default CSV parameters.
    sep:
        Field delimiter of the CSV files. Defaults to ``cfg.csv_sep``.
    encoding:
        Text encoding of the CSV files. Defaults to ``cfg.csv_encoding``.

    """
    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding
    try:
        df = pd.read_csv(
            input_path,
            sep=sep,
            encoding=encoding,
            dtype=str,
            dtype_backend="pyarrow",
        )
    except ImportError:  # pragma: no cover - pyarrow optional
        df = pd.read_csv(input_path, sep=sep, encoding=encoding, dtype=str)
    processed = postprocess_assays(df)
    processed.to_csv(output_path, index=False, sep=sep, encoding=encoding)
