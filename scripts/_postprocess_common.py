"""Shared helpers for post-processing CLI entry points."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd

from library import io
from library.common.csv_utils import write_csv_chunks_deterministic
from library.postprocessing.common.config import PipelineConfig, load_pipeline_config
from library.postprocessing.common.logging import (
    PipelineRunMetrics,
    build_report_payload,
    dump_report,
)
from library.postprocessing.common.schema import DataFrameSchema
from library.postprocessing.common.types import SchemaValidationError
from library.postprocessing.common.utils import collect_postprocess_metrics

LOG_DIR_ENV = "CHEMBL_POSTPROCESS_LOG_DIR"
DEFAULT_LOG_DIR = Path("logs")


@dataclass(slots=True)
class CsvRuntimeConfig:
    """Runtime configuration for CSV I/O operations."""

    sep: str
    encoding: str
    chunksize: int = 10000


def event_prefix(table: str) -> str:
    """Return the structured logging prefix for ``table``."""

    return f"{table}_postprocess"


def get_pipeline_config(table: str, config_path: Path | None) -> PipelineConfig:
    """Load the pipeline configuration for ``table``."""

    return load_pipeline_config(table, config_path)


def get_csv_runtime_config(config: PipelineConfig) -> CsvRuntimeConfig:
    """Derive CSV runtime parameters from ``config``."""

    params = dict(config.params or {})
    io_params = dict(params.get("io", {}))
    sep = str(io_params.get("csv_sep", ","))
    encoding = str(io_params.get("encoding", "utf-8-sig"))
    chunksize = int(io_params.get("chunksize", 10000))
    return CsvRuntimeConfig(sep=sep, encoding=encoding, chunksize=chunksize)


def get_default_log_level(config: PipelineConfig) -> str:
    """Return the default log level declared in ``config``."""

    params = dict(config.params or {})
    defaults = params.get("defaults", {}) or {}
    level = defaults.get("log_level", "INFO")
    return str(level).upper()


def ensure_pipeline_version_column(
    df: pd.DataFrame, pipeline_version: str | None
) -> pd.DataFrame:
    """Ensure ``pipeline_version`` is represented as a string column."""

    if not pipeline_version:
        return df

    prepared = df.copy(deep=True)
    version = str(pipeline_version)
    if "pipeline_version" in prepared.columns:
        series = prepared["pipeline_version"].astype("string")
        prepared["pipeline_version"] = series.fillna(version)
    else:
        prepared["pipeline_version"] = pd.Series(
            pd.array([version] * len(prepared.index), dtype="string"),
            index=prepared.index,
        )
    return prepared


def load_input_frame(
    table: str,
    path: Path,
    csv_cfg: CsvRuntimeConfig,
    *,
    logger,
) -> pd.DataFrame:
    """Load the raw output frame for ``table`` from ``path``."""

    prefix = event_prefix(table)
    logger.info(
        f"{prefix}_load_start",
        input=str(path),
        encoding=csv_cfg.encoding,
        separator=csv_cfg.sep,
    )
    namespace = SimpleNamespace(csv_sep=csv_cfg.sep, csv_encoding=csv_cfg.encoding)
    frame = io.read_csv(path, cfg=namespace)
    rows, cols = frame.shape
    logger.info(
        f"{prefix}_load_done",
        rows=int(rows),
        columns=int(cols),
    )
    return frame


def run_postprocess_steps(
    table: str,
    df: pd.DataFrame,
    runner: Callable[..., tuple[pd.DataFrame, PipelineRunMetrics]],
    pipeline_version: str | None,
    *,
    logger,
) -> tuple[pd.DataFrame, PipelineRunMetrics]:
    """Execute postprocessing ``runner`` for ``table``."""

    prefix = event_prefix(table)
    logger.info(
        f"{prefix}_pipeline_start",
        rows=int(df.shape[0]),
        columns=int(df.shape[1]),
        pipeline_version=pipeline_version,
    )
    processed, metrics = runner(
        df,
        pipeline_version=pipeline_version,
        logger=logger,
    )
    summary = metrics.summary()
    logger.info(f"{prefix}_pipeline_done", **summary)
    return processed, metrics


def validate_postprocess_frame(
    table: str,
    df: pd.DataFrame,
    validator: Callable[..., pd.DataFrame],
    schema: DataFrameSchema,
    pipeline_version: str | None,
    *,
    logger,
) -> pd.DataFrame:
    """Validate ``df`` using ``validator`` aligned with ``schema``."""

    prefix = event_prefix(table)
    logger.info(
        f"{prefix}_validate_start",
        schema=schema.__class__.__name__,
    )
    prepared = ensure_pipeline_version_column(df, pipeline_version)
    validated = validator(prepared, context=f"{table}_postprocess")
    logger.info(
        f"{prefix}_validate_done",
        rows=int(validated.shape[0]),
        columns=int(validated.shape[1]),
    )
    return validated


def _resolve_key_columns(schema: DataFrameSchema, df: pd.DataFrame) -> Sequence[str]:
    candidates: Sequence[str] = schema.sort_by or schema.required_columns or ()
    selected = [column for column in candidates if column in df.columns]
    if not selected:
        remaining = [column for column in df.columns if column]
        if not remaining:
            raise SchemaValidationError(
                "export",
                "no columns available to determine deterministic ordering",
            )
        selected = [remaining[0]]
    return tuple(selected)


def export_postprocess_frame(
    table: str,
    df: pd.DataFrame,
    output_path: Path,
    csv_cfg: CsvRuntimeConfig,
    schema: DataFrameSchema,
    *,
    logger,
) -> Path:
    """Persist ``df`` deterministically to ``output_path``."""

    prefix = event_prefix(table)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(
        f"{prefix}_export_start",
        output=str(output_path),
        rows=int(df.shape[0]),
        columns=int(df.shape[1]),
        encoding=csv_cfg.encoding,
        separator=csv_cfg.sep,
    )

    key_columns = _resolve_key_columns(schema, df)
    column_order: Sequence[str] | None
    if schema.column_order:
        column_order = tuple(col for col in schema.column_order if col in df.columns)
    else:
        column_order = tuple(df.columns)

    write_csv_chunks_deterministic(
        (df,),
        output_path,
        col_order=column_order,
        key_cols=key_columns,
        sep=csv_cfg.sep,
        encoding=csv_cfg.encoding,
        chunksize=csv_cfg.chunksize,
    )

    logger.info(f"{prefix}_export_done", output=str(output_path))
    return output_path


def generate_metrics_report(
    table: str,
    output_path: Path,
    csv_cfg: CsvRuntimeConfig,
    runner: Callable[..., tuple[pd.DataFrame, PipelineRunMetrics]],
    *,
    pipeline_version: str | None,
    extras: Mapping[str, Any] | None,
    logger,
    pipeline_metrics: PipelineRunMetrics | None = None,
) -> tuple[PipelineRunMetrics | None, Path | None]:
    """Create a structured metrics report for ``table``."""

    prefix = event_prefix(table)
    if pipeline_metrics is not None:
        report_path = output_path.parent / f"{table}.postprocess.report.json"
        payload = build_report_payload(
            table=table,
            metrics=pipeline_metrics,
            output_path=str(output_path),
            extras=extras,
        )
        dump_report(report_path, payload)
        logger.info(
            f"{prefix}_report_written",
            report=str(report_path),
            rows=pipeline_metrics.output_rows,
            columns=pipeline_metrics.output_columns,
        )
        return pipeline_metrics, report_path

    metrics, report_path = collect_postprocess_metrics(
        table=table,
        output_path=output_path,
        csv_sep=csv_cfg.sep,
        csv_encoding=csv_cfg.encoding,
        output_dir=output_path.parent,
        runner=runner,
        logger=logger,
        pipeline_version=pipeline_version,
        report_extras=extras or {},
    )
    if metrics is not None and report_path is not None:
        logger.info(
            f"{prefix}_report_written",
            report=str(report_path),
            rows=metrics.output_rows,
            columns=metrics.output_columns,
        )
    else:
        logger.info(f"{prefix}_report_skipped", output=str(output_path))
    return metrics, report_path


__all__ = [
    "CsvRuntimeConfig",
    "DEFAULT_LOG_DIR",
    "LOG_DIR_ENV",
    "ensure_pipeline_version_column",
    "event_prefix",
    "export_postprocess_frame",
    "generate_metrics_report",
    "get_csv_runtime_config",
    "get_default_log_level",
    "get_pipeline_config",
    "load_input_frame",
    "run_postprocess_steps",
    "validate_postprocess_frame",
]
