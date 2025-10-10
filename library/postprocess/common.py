from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Callable, Iterable, Mapping, Any

import pandas as pd

from library import io
from library.common.csv_utils import write_csv_chunks_deterministic
from library.postprocessing.activities import ACTIVITY_SCHEMA, run_activity_pipeline
from library.postprocessing.assays import ASSAY_SCHEMA, run_assay_pipeline
from library.postprocessing.common.config import PipelineConfig, load_pipeline_config
from library.postprocessing.common.logging import (
    PipelineRunMetrics,
    build_report_payload,
    dump_report,
)
from library.postprocessing.common.schema import DataFrameSchema
from library.postprocessing.common.types import SchemaValidationError
from library.postprocessing.common.utils import collect_postprocess_metrics
from library.postprocessing.documents import DOCUMENT_SCHEMA, run_document_pipeline
from library.postprocessing.targets import TARGET_SCHEMA, run_target_pipeline
from library.postprocessing.testitem import TESTITEM_SCHEMA, run_testitem_pipeline

__all__ = ["PostprocessResult", "SUPPORTED_TABLES", "run_postprocessing_pipeline"]


_LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class _DomainResources:
    """Bundle runtime artefacts required to run a post-processing pipeline."""

    runner: Callable[[pd.DataFrame], tuple[pd.DataFrame, PipelineRunMetrics]]
    schema: DataFrameSchema


_DOMAIN_RESOURCES: dict[str, _DomainResources] = {
    "activities": _DomainResources(run_activity_pipeline, ACTIVITY_SCHEMA),
    "assays": _DomainResources(run_assay_pipeline, ASSAY_SCHEMA),
    "documents": _DomainResources(run_document_pipeline, DOCUMENT_SCHEMA),
    "targets": _DomainResources(run_target_pipeline, TARGET_SCHEMA),
    "testitems": _DomainResources(run_testitem_pipeline, TESTITEM_SCHEMA),
}

SUPPORTED_TABLES: tuple[str, ...] = tuple(sorted(_DOMAIN_RESOURCES))

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LOG_DIR: Path = (_PROJECT_ROOT / "logs").resolve()
LOG_DIR_ENV = "CHEMBL_POSTPROCESS_LOG_DIR"


@dataclass(frozen=True)
class CsvRuntimeConfig:
    """Runtime CSV parameters derived from the domain configuration."""

    separator: str
    encoding: str
    chunksize: int


@dataclass(frozen=True)
class PostprocessResult:
    """Container returned by :func:`run_postprocessing_pipeline`."""

    table: str
    output_path: Path
    report_path: Path
    metrics: PipelineRunMetrics


def _event_prefix(table: str) -> str:
    return f"{table}_postprocess"


def _resolve_domain(table: str) -> _DomainResources:
    try:
        return _DOMAIN_RESOURCES[table]
    except KeyError as exc:  # pragma: no cover - defensive guard
        raise ValueError(f"unsupported postprocess table: {table!r}") from exc


def _load_input_frame(table: str, path: Path, csv_cfg: CsvRuntimeConfig, logger: logging.Logger) -> pd.DataFrame:
    prefix = _event_prefix(table)
    logger.info(
        f"{prefix}_load_start",
        input=str(path),
        encoding=csv_cfg.encoding,
        separator=csv_cfg.separator,
    )
    namespace = SimpleNamespace(csv_sep=csv_cfg.separator, csv_encoding=csv_cfg.encoding)
    frame = io.read_csv(path, cfg=namespace)
    rows, cols = frame.shape
    logger.info(f"{prefix}_load_done", rows=int(rows), columns=int(cols))
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

    prefix = _event_prefix(table)
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

    prefix = _event_prefix(table)
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


def _resolve_key_columns(schema: DataFrameSchema, df: pd.DataFrame) -> tuple[str, ...]:
    candidates = schema.sort_by or schema.required_columns or ()
    selected = [column for column in candidates if column in df.columns and column]
    if selected:
        return tuple(selected)
    remaining = [column for column in df.columns if column]
    if not remaining:
        raise SchemaValidationError(
            "postprocess",
            "no columns available to determine ordering",
        )
    return (remaining[0],)


def _write_output_frame(
    table: str,
    df: pd.DataFrame,
    output_path: Path,
    csv_cfg: CsvRuntimeConfig,
    schema: DataFrameSchema,
    logger: logging.Logger,
) -> Path:
    """Persist ``df`` deterministically to ``output_path``."""

    prefix = _event_prefix(table)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(
        f"{prefix}_export_start",
        output=str(output_path),
        rows=int(df.shape[0]),
        columns=int(df.shape[1]),
        encoding=csv_cfg.encoding,
        separator=csv_cfg.separator,
    )
    if schema.column_order:
        column_order = tuple(col for col in schema.column_order if col in df.columns)
    else:
        column_order = tuple(df.columns)
    key_columns = _resolve_key_columns(schema, df)
    write_csv_chunks_deterministic(
        (df,),
        output_path,
        col_order=column_order,
        key_cols=key_columns,
        sep=csv_cfg.separator,
        encoding=csv_cfg.encoding,
        chunksize=csv_cfg.chunksize,
    )
    logger.info(f"{prefix}_export_done", output=str(output_path))
    return output_path


@dataclass(slots=True)
class PostprocessingPipelineConfig:
    """Configuration bundle required to execute a postprocessing pipeline."""

    pipeline_config: "PipelineConfig"
    csv_runtime_config: "CsvRuntimeConfig"
    runner: Callable[..., tuple["pd.DataFrame", "PipelineRunMetrics"]]
    validator: Callable[..., "pd.DataFrame"]
    schema: "DataFrameSchema"
    logger: Any
    pipeline_version: str | None = None

@dataclass(slots=True)
class PostprocessingPipelineResult:
    """Result payload returned by :func:`run_postprocessing_pipeline`."""

    dataframe: "pd.DataFrame"
    metrics: "PipelineRunMetrics | None"
    output_path: Path
    report_path: Path | None

def _write_metrics_report(
    table: str,
    metrics: "PipelineRunMetrics",
    output_path: Path,
    logger: logging.Logger,
) -> Path:
    report_path = output_path.parent / f"{table}.postprocess.report.json"
    payload = build_report_payload(
        table=table,
        metrics=metrics,
        output_path=str(output_path),
    )
    dump_report(report_path, payload)
    logger.info(f"{_event_prefix(table)}_report_written", report=str(report_path))
    return report_path

def run_postprocessing_pipeline(
    table_name: str,
    input_path: Path | str,
    output_path: Path | str,
    config: PostprocessingPipelineConfig,
) -> PostprocessingPipelineResult:
    """Execute the canonical load → process → validate → save pipeline."""

    csv_cfg = config.csv_runtime_config
    logger = config.logger
    pipeline_config = config.pipeline_config
    runner = config.runner
    validator = config.validator
    schema = config.schema

    resolved_input = Path(input_path)
    resolved_output = Path(output_path)

    if not resolved_input.exists():
        raise FileNotFoundError(resolved_input)

    pipeline_version = config.pipeline_version or pipeline_config.pipeline_version

    frame = load_input_frame(table_name, resolved_input, csv_cfg, logger=logger)
    processed, metrics = run_postprocess_steps(
        table_name,
        frame,
        runner,
        pipeline_version,
        logger=logger,
    )

    effective_version = metrics.pipeline_version if metrics else pipeline_version
    validated = validate_postprocess_frame(
        table_name,
        processed,
        validator,
        schema,
        effective_version,
        logger=logger,
    )
    export_postprocess_frame(
        table_name,
        validated,
        resolved_output,
        csv_cfg,
        schema,
        logger=logger,
    )

    report_path: Path | None = None

    # Write metrics report after exporting CSV
    if metrics is not None:
        report_path = _write_metrics_report(
            table_name, metrics, resolved_output, logger
        )

    return PostprocessingPipelineResult(
        dataframe=validated,
        metrics=metrics,
        output_path=resolved_output,
        report_path=report_path,
    )

def generate_metrics_report(
    table: str,
    output_path: Path,
    csv_cfg: "CsvRuntimeConfig",
    runner: Callable[..., tuple["pd.DataFrame", "PipelineRunMetrics"]],
    *,
    pipeline_version: str | None,
    extras: Mapping[str, Any] | None,
    logger,
    pipeline_metrics: "PipelineRunMetrics | None" = None,
) -> tuple["PipelineRunMetrics | None", Path | None]:
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
        csv_sep=csv_cfg.separator,
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
    "PostprocessingPipelineConfig",
    "PostprocessingPipelineResult",
    "ensure_pipeline_version_column",
    "event_prefix",
    "export_postprocess_frame",
    "generate_metrics_report",
    "get_csv_runtime_config",
    "get_default_log_level",
    "get_pipeline_config",
    "load_input_frame",
    "run_postprocess_steps",
    "run_postprocessing_pipeline",
    "validate_postprocess_frame",
]


def event_prefix(table: str) -> str:
    """Return the structured logging prefix for ``table`` postprocess events."""

    return _event_prefix(table)


def ensure_pipeline_version_column(
    df: pd.DataFrame, pipeline_version: str | None
) -> pd.DataFrame:
    """Ensure that ``df`` exposes a string ``pipeline_version`` column."""

    prepared = df.copy(deep=True)

    if "pipeline_version" in prepared.columns:
        prepared["pipeline_version"] = prepared["pipeline_version"].astype("string")
    else:
        prepared.insert(
            len(prepared.columns),
            "pipeline_version",
            pd.Series(pd.NA, index=prepared.index, dtype="string"),
        )

    if pipeline_version:
        prepared["pipeline_version"] = (
            prepared["pipeline_version"].fillna(str(pipeline_version)).astype("string")
        )

    return prepared


def load_input_frame(
    table: str, path: Path | str, csv_cfg: CsvRuntimeConfig, *, logger: logging.Logger
) -> pd.DataFrame:
    """Read ``path`` into a DataFrame honouring ``csv_cfg`` options."""

    return _load_input_frame(table, Path(path), csv_cfg, logger)


def export_postprocess_frame(
    table: str,
    df: pd.DataFrame,
    output_path: Path | str,
    csv_cfg: CsvRuntimeConfig,
    schema: DataFrameSchema,
    *,
    logger: logging.Logger,
) -> Path:
    """Persist ``df`` deterministically to ``output_path``."""

    return _write_output_frame(table, df, Path(output_path), csv_cfg, schema, logger)


def get_pipeline_config(table: str, override: Path | str | None = None) -> PipelineConfig:
    """Return the declarative pipeline configuration for ``table``."""

    resolved_override = None if override is None else Path(override)
    return load_pipeline_config(table, resolved_override)


def get_csv_runtime_config(pipeline_config: PipelineConfig) -> CsvRuntimeConfig:
    """Derive runtime CSV parameters from ``pipeline_config``."""

    params = dict(pipeline_config.params or {})
    io_params = dict(params.get("io", {}) or {})
    separator = str(io_params.get("csv_sep", ","))
    encoding = str(io_params.get("encoding", "utf-8-sig"))
    chunksize = int(io_params.get("chunksize", 10000))
    return CsvRuntimeConfig(separator=separator, encoding=encoding, chunksize=chunksize)


def get_default_log_level(pipeline_config: PipelineConfig) -> str:
    """Return the default log level declared in ``pipeline_config``."""

    params = dict(pipeline_config.params or {})
    defaults = dict(params.get("defaults", {}) or {})
    level = defaults.get("log_level")
    if isinstance(level, str) and level.strip():
        return level.strip()
    return "INFO"
