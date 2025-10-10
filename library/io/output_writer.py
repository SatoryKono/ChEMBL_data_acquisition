"""High-level helpers to persist standard pipeline artefacts deterministically."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from collections.abc import Sequence

import pandas as pd

from ..common.csv_utils import write_csv_deterministic
from ..config import IoCfg


@dataclass(frozen=True, slots=True)
class StandardOutputArtifacts:
    """Paths of the primary dataset and accompanying QA reports."""

    dataset: Path
    quality_report: Path
    correlation_report: Path


def _ensure_output_directory(path: Path, *, cfg: IoCfg) -> None:
    """Create ``path`` when missing if permitted by :class:`IoCfg`."""

    if path.exists():
        if not path.is_dir():
            raise NotADirectoryError(f"{path} is not a directory")
        return
    if cfg.exist_ok:
        path.mkdir(parents=True, exist_ok=True)
    else:
        raise FileNotFoundError(f"{path} does not exist")


def save_standard_outputs(
    dataset: pd.DataFrame,
    quality_report: pd.DataFrame,
    correlation_report: pd.DataFrame,
    *,
    table_name: str,
    date_tag: str,
    cfg: IoCfg,
    key_columns: Sequence[str] | None = None,
) -> StandardOutputArtifacts:
    """Persist canonical output artefacts and return their locations.

    Parameters
    ----------
    dataset:
        Final merged dataset for the pipeline.
    quality_report:
        Quality control summary table.
    correlation_report:
        Data correlation insights produced during QA.
    table_name:
        Resolved name of the output table (e.g. ``"document"``).
    date_tag:
        Deterministic UTC date string appended to filenames.
    cfg:
        Active IO configuration controlling output directory handling.
    key_columns:
        Optional sequence of columns used for deterministic ordering of the
        dataset export. When omitted the first column of ``dataset`` is used.

    Returns
    -------
    StandardOutputArtifacts
        Absolute paths for the dataset, quality report and correlation report.
    """

    output_dir = Path(cfg.output_dir)
    _ensure_output_directory(output_dir, cfg=cfg)

    stem = f"output.{table_name}_{date_tag}"
    dataset_path = output_dir / f"{stem}.csv"
    quality_path = output_dir / f"{stem}_quality_report_table.csv"
    correlation_path = output_dir / f"{stem}_data_correlation_report_table.csv"

    key_cols = list(key_columns or [])
    if not key_cols and not dataset.empty:
        key_cols = [str(dataset.columns[0])]

    write_csv_deterministic(
        dataset,
        dataset_path,
        key_cols=key_cols,
        cfg=cfg,
    )
    write_csv_deterministic(
        quality_report,
        quality_path,
        key_cols=list(quality_report.columns[:1]),
        cfg=cfg,
    )
    corr_key_cols: list[str] = []
    if not correlation_report.empty:
        corr_key_cols = [str(correlation_report.columns[0])]
    write_csv_deterministic(
        correlation_report,
        correlation_path,
        key_cols=corr_key_cols,
        cfg=cfg,
    )

    return StandardOutputArtifacts(
        dataset=dataset_path,
        quality_report=quality_path,
        correlation_report=correlation_path,
    )
