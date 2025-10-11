"""High-level helpers to persist standard pipeline artefacts deterministically."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import pandas as pd

from ..common.csv_utils import write_csv_deterministic
from ..config import IoCfg

__all__ = ["StandardOutputArtifacts", "save_standard_outputs"]


@dataclass(frozen=True, slots=True)
class StandardOutputArtifacts:
    """Paths of the primary dataset and accompanying QA reports."""

    dataset: Path
    correlation_report: Path
    quality_report: Path


def _ensure_output_directory(path: Path, *, exist_ok: bool) -> None:
    """Ensure ``path`` exists honouring ``exist_ok`` semantics."""

    if path.exists():
        if not path.is_dir():  # pragma: no cover - defensive guard
            raise NotADirectoryError(f"{path} is not a directory")
        if not exist_ok:
            raise FileExistsError(f"{path} already exists")
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
    io_cfg: IoCfg | None = None,
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
    io_cfg:
        Optional :class:`~library.config.IoCfg` controlling CSV formatting and
        destination directory. When omitted a fresh instance with default
        values is used, honouring environment overrides such as
        ``CHEMBL_DA_BASE_PATH``.

    Returns
    -------
    StandardOutputArtifacts
        Paths to the dataset, correlation report and quality report.
    """

    cfg = io_cfg or IoCfg()
    output_dir = Path(cfg.output_dir)

    _ensure_output_directory(output_dir, exist_ok=cfg.exist_ok)

    stem = f"output.{table_name}_{date_tag}"
    dataset_path = output_dir / f"{stem}.csv"
    correlation_path = output_dir / f"{stem}_data_correlation_report_table.csv"
    quality_path = output_dir / f"{stem}_quality_report_table.csv"

    key_cols = list(key_columns or [])
    if not key_cols and not df_main.empty:
        key_cols = [str(df_main.columns[0])]

    write_csv_deterministic(
        df_main,
        dataset_path,
        key_cols=key_cols,
        sep=cfg.csv_sep,
        encoding=cfg.csv_encoding,
        chunksize=cfg.csv_chunksize,
        sort_chunksize=cfg.csv_chunksize,
    )

    qc_key_cols: list[str] = []
    if not df_qc.empty:
        qc_key_cols = [str(df_qc.columns[0])]
    write_csv_deterministic(
        df_qc,
        quality_path,
        key_cols=qc_key_cols,
        sep=cfg.csv_sep,
        encoding=cfg.csv_encoding,
        chunksize=cfg.csv_chunksize,
        sort_chunksize=cfg.csv_chunksize,
    )

    corr_key_cols: list[str] = []
    if not df_corr.empty:
        corr_key_cols = [str(df_corr.columns[0])]
    write_csv_deterministic(
        df_corr,
        correlation_path,
        key_cols=corr_key_cols,
        sep=cfg.csv_sep,
        encoding=cfg.csv_encoding,
        chunksize=cfg.csv_chunksize,
        sort_chunksize=cfg.csv_chunksize,
    )

    return StandardOutputArtifacts(
        dataset=dataset_path,
        correlation_report=correlation_path,
        quality_report=quality_path,
    )
