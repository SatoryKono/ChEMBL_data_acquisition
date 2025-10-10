"""Pipeline orchestration utilities for postprocessing transformations."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import pandas as pd

from . import logging as logging_utils


def infer_pipeline_version(frame: pd.DataFrame) -> str | None:
    """Return the first non-empty ``pipeline_version`` value from ``frame``."""

    if "pipeline_version" not in frame.columns or frame.empty:
        return None
    series = frame["pipeline_version"].dropna()
    if series.empty:
        return None
    value = str(series.iloc[0]).strip()
    return value or None


def collect_postprocess_metrics(
    *,
    table: str,
    output_path: Path,
    csv_sep: str | None,
    csv_encoding: str | None,
    output_dir: Path | str,
    runner: Callable[..., tuple[pd.DataFrame, logging_utils.PipelineRunMetrics]],
    logger,
    pipeline_version: str | None = None,
    report_extras: Mapping[str, Any] | None = None,
) -> tuple[logging_utils.PipelineRunMetrics | None, Path | None]:
    """Load ``output_path`` and generate a postprocess metrics report."""

    event_prefix = f"{table}_postprocess"
    if not output_path.exists():
        logger.warning(
            f"{event_prefix}_report_missing_output",
            output_postprocessed=str(output_path),
        )
        return None, None

    read_kwargs: dict[str, Any] = {}
    if csv_sep is not None:
        read_kwargs["sep"] = csv_sep
    if csv_encoding is not None:
        read_kwargs["encoding"] = csv_encoding

    try:
        frame = pd.read_csv(output_path, **read_kwargs)
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning(
            f"{event_prefix}_report_load_failed",
            output_postprocessed=str(output_path),
            error=str(exc),
        )
        return None, None

    effective_version = pipeline_version or infer_pipeline_version(frame)

    try:
        _, metrics = runner(
            frame,
            pipeline_version=effective_version,
            logger=logger,
        )
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning(
            f"{event_prefix}_report_generation_failed",
            output_postprocessed=str(output_path),
            error=str(exc),
        )
        return None, None

    report_dir = Path(output_dir)
    report_path = report_dir / f"{table}.postprocess.report.json"
    payload = logging_utils.build_report_payload(
        table=table,
        metrics=metrics,
        output_postprocessed=str(output_path),
        extras=report_extras,
    )
    logging_utils.dump_report(report_path, payload)
    return metrics, report_path


__all__ = ["collect_postprocess_metrics", "infer_pipeline_version"]
