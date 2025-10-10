from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Callable, Iterable

import pandas as pd

from library import io
from library.common.csv_utils import write_csv_chunks_deterministic
from library.postprocessing.activities import ACTIVITY_SCHEMA, run_activity_pipeline
from library.postprocessing.assays import ASSAY_SCHEMA, run_assay_pipeline
from library.postprocessing.common.config import load_pipeline_config
from library.postprocessing.common.logging import (
    PipelineRunMetrics,
    build_report_payload,
    dump_report,
)
from library.postprocessing.common.schema import DataFrameSchema
from library.postprocessing.common.types import SchemaValidationError
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


def _load_csv_runtime_config(table: str, override: Path | None) -> tuple[CsvRuntimeConfig, str | None]:
    config = load_pipeline_config(table, override)
    params = dict(config.params or {})
    io_params = dict(params.get("io", {}) or {})
    separator = str(io_params.get("csv_sep", ","))
    encoding = str(io_params.get("encoding", "utf-8-sig"))
    chunksize = int(io_params.get("chunksize", 10000))
    return CsvRuntimeConfig(separator=separator, encoding=encoding, chunksize=chunksize), config.pipeline_version


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


def _resolve_key_columns(schema: DataFrameSchema, df: pd.DataFrame) -> Iterable[str]:
    candidates = schema.sort_by or schema.required_columns or ()
    selected = [column for column in candidates if column in df.columns and column]
    if selected:
        return tuple(selected)
    remaining = [column for column in df.columns if column]
    if not remaining:
        raise SchemaValidationError("postprocess", "no columns available to determine ordering")
    return (remaining[0],)


def _write_output_frame(
    table: str,
    df: pd.DataFrame,
    output_path: Path,
    csv_cfg: CsvRuntimeConfig,
    schema: DataFrameSchema,
    logger: logging.Logger,
) -> Path:
    prefix = _event_prefix(table)
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


def _write_metrics_report(
    table: str,
    metrics: PipelineRunMetrics,
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
    table: str,
    input_path: Path | str,
    *,
    output_path: Path | str | None = None,
    config_path: Path | str | None = None,
    logger: logging.Logger | None = None,
) -> PostprocessResult:
    """Run the post-processing pipeline for ``table`` on ``input_path``.

    Parameters
    ----------
    table:
        Domain identifier (for example ``"assays"`` or ``"targets"``).
    input_path:
        Path to the freshly generated CSV artefact.
    output_path:
        Optional override for the post-processed CSV destination. When omitted
        the file is written next to ``input_path`` under the
        ``output_postprocessed.<table>.csv`` name.
    config_path:
        Optional override for the YAML pipeline configuration. Defaults to the
        canonical ``config/pipeline/<table>.yaml`` location.
    logger:
        Logger used for structured progress messages. The module logger is used
        when omitted.

    Returns
    -------
    PostprocessResult
        Object containing the destination CSV path, the generated metrics
        report location and the collected pipeline metrics.
    """

    resources = _resolve_domain(table)
    log = logger or _LOGGER
    input_path = Path(input_path)
    if output_path is None:
        output_path = input_path.with_name(f"output_postprocessed.{table}.csv")
    else:
        output_path = Path(output_path)
    csv_cfg, pipeline_version_override = _load_csv_runtime_config(
        table, Path(config_path) if config_path is not None else None
    )
    frame = _load_input_frame(table, input_path, csv_cfg, log)
    processed, metrics = resources.runner(
        frame,
        pipeline_version=pipeline_version_override,
        logger=log,
    )
    output_path = _write_output_frame(
        table,
        processed,
        output_path,
        csv_cfg,
        resources.schema,
        log,
    )
    report_path = _write_metrics_report(table, metrics, output_path, log)
    log.info(
        f"{_event_prefix(table)}_done",
        output=str(output_path),
        report=str(report_path),
        rows=metrics.output_rows,
        columns=metrics.output_columns,
    )
    return PostprocessResult(
        table=table,
        output_path=output_path,
        report_path=report_path,
        metrics=metrics,
    )

