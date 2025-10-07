"""Assay table post-processing utilities.

This module provides helper functions to transform assay data retrieved from
ChEMBL.  The main feature adds an ``assay_with_same_target`` column that counts
how many assays share the same ``document_chembl_id`` and
``target_chembl_id`` combination.  The transformation mimics the behaviour of a
Power Query script used in the original workflow.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.schemas import AssayPostprocessSchema

from ...common.csv_utils import write_csv_deterministic
from ...common.log import logger
from ...common.pandas_utils import read_csv_pyarrow
from ...config import IoCfg


def postprocess_assays(df: pd.DataFrame) -> pd.DataFrame:
    """Augment assay data with per-target assay counts.

    The function adds a new column named ``assay_with_same_target`` to ``df``.  The
    value represents how many assays share the same ``document_chembl_id`` and
    ``target_chembl_id`` pair as the current row.

    Parameters
    ----------
    df:
        DataFrame produced by :func:`library.integration.chembl_library.get_assays`.

    Returns
    -------
    pandas.DataFrame
        Copy of the input with the additional ``assay_with_same_target`` column.

    """
    df = df.copy()
    required_columns = ["document_chembl_id", "target_chembl_id"]

    # ``pandera`` raises ``SchemaError`` on completely empty frames because the
    # required columns are absent.  When the pipeline processes an empty chunk we
    # still want to synthesise the expected columns so that downstream hooks see
    # a consistent schema.  Additionally, normalise the dtype to ``string`` so
    # that subsequent concatenation does not introduce ``object`` columns.
    for column in required_columns:
        if column not in df.columns:
            df[column] = pd.Series(dtype="string")
        else:
            df[column] = df[column].astype("string")

    if df.empty:
        ordered = required_columns + [
            column for column in df.columns if column not in required_columns
        ]
        df = df.reindex(columns=ordered)
    else:
        df = AssayPostprocessSchema.validate(df)
    result = df.copy()

    if result.empty:
        result["assay_with_same_target"] = pd.Series(dtype="Int64")
        return result

    group_cols = ["document_chembl_id", "target_chembl_id"]
    groups = result.groupby(group_cols)
    logger.debug("Calculated counts for %d document/target groups", groups.ngroups)
    counts = groups["document_chembl_id"].transform("size").astype("Int64")
    result["assay_with_same_target"] = counts
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
        Path to the CSV file produced by ``scripts/get_assay_data.py``.
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
        df = read_csv_pyarrow(
            input_path,
            sep=sep,
            encoding=encoding,
            dtype=str,
        )
    except ImportError:  # pragma: no cover - pyarrow optional
        df = pd.read_csv(input_path, sep=sep, encoding=encoding, dtype=str)
    processed = postprocess_assays(df)
    write_csv_deterministic(
        processed,
        output_path,
        col_order=list(processed.columns),
        key_cols=["assay_chembl_id"],
        sep=sep,
        encoding=encoding,
        cfg=None,
    )
