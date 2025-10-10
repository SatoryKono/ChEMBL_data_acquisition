"""High-level helpers to persist standard pipeline artefacts deterministically."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

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

    def _write(frame: pd.DataFrame, path: Path) -> None:
        key_cols = sorted(frame.columns)
        col_order = list(frame.columns)
        write_csv_deterministic(
            frame,
            path,
            key_cols=key_cols,
            col_order=col_order if col_order else None,
            sep=cfg.csv_sep,
            encoding=cfg.csv_encoding,
            cfg=cfg,
        )

    _write(dataset, dataset_path)
    _write(quality_report, quality_path)
    _write(correlation_report, correlation_path)

    return StandardOutputArtifacts(
        dataset=dataset_path,
        quality_report=quality_path,
        correlation_report=correlation_path,
    )
