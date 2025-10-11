"""High-level helpers to persist standard pipeline artefacts deterministically."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import pandas as pd

from ..common.csv_utils import write_csv_deterministic


OUTPUT_DIR = Path("data/output")


@dataclass(frozen=True, slots=True)
class StandardOutputArtifacts:
    """Paths of the primary dataset and accompanying QA reports."""

    dataset: Path
    correlation_report: Path
    quality_report: Path


def _ensure_output_directory(path: Path) -> None:
    """Create ``path`` when missing."""

    if path.exists():
        if not path.is_dir():  # pragma: no cover - defensive guard
            raise NotADirectoryError(f"{path} is not a directory")
        return
    path.mkdir(parents=True, exist_ok=True)


def save_standard_outputs(
    df_main: pd.DataFrame,
    df_corr: pd.DataFrame,
    df_qc: pd.DataFrame,
    table_name: str,
    date_tag: str,
    *,
    key_columns: Sequence[str] | None = None,
    dataset_path: Path | None = None,
    output_dir: Path | None = None,
) -> StandardOutputArtifacts:
    """Persist the canonical dataset together with QC artefacts.

    Parameters
    ----------
    df_main:
        Aggregated dataset for the table.
    df_corr:
        Correlation metrics generated for ``df_main``.
    df_qc:
        Quality-control summary for ``df_main``.
    table_name:
        Logical name of the table (``document``, ``target`` and so on).
    date_tag:
        UTC date token formatted as ``YYYYMMDD``.
    key_columns:
        Optional sequence of columns used to deterministically order
        ``df_main`` before writing.
    dataset_path:
        Explicit path for the dataset CSV. When provided the correlation and
        quality reports are emitted alongside it using the same stem.
    output_dir:
        Directory where artefacts should be written. Ignored when
        ``dataset_path`` is provided. Defaults to ``data/output``.

    Returns
    -------
    StandardOutputArtifacts
        Paths to the dataset, correlation report and quality report.
    """

    if dataset_path is not None and output_dir is not None:
        output_dir = Path(output_dir)
        dataset_path = Path(dataset_path)
        if dataset_path.parent != output_dir:
            raise ValueError(
                "dataset_path parent must match output_dir when both are provided"
            )

    if dataset_path is not None:
        dataset_path = Path(dataset_path)
        base_dir = dataset_path.parent
        stem = dataset_path.stem
    else:
        base_dir = Path(output_dir) if output_dir is not None else OUTPUT_DIR
        stem = f"output.{table_name}_{date_tag}"
        dataset_path = base_dir / f"{stem}.csv"

    _ensure_output_directory(base_dir)

    correlation_path = base_dir / f"{stem}_data_correlation_report_table.csv"
    quality_path = base_dir / f"{stem}_quality_report_table.csv"

    key_cols = list(key_columns or [])
    if not key_cols and not df_main.empty:
        key_cols = [str(df_main.columns[0])]

    write_csv_deterministic(df_main, dataset_path, key_cols=key_cols)

    qc_key_cols: list[str] = []
    if not df_qc.empty:
        qc_key_cols = [str(df_qc.columns[0])]
    write_csv_deterministic(df_qc, quality_path, key_cols=qc_key_cols)

    corr_key_cols: list[str] = []
    if not df_corr.empty:
        corr_key_cols = [str(df_corr.columns[0])]
    write_csv_deterministic(df_corr, correlation_path, key_cols=corr_key_cols)

    return StandardOutputArtifacts(
        dataset=dataset_path,
        correlation_report=correlation_path,
        quality_report=quality_path,
    )
