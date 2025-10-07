"""Pipeline orchestration utilities for postprocessing transformations."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Any, Callable, Mapping
import pandas as pd

from . import logging as logging_utils
from .io import clone_dataframe, ensure_dataframe
from .schema import DataFrameSchema, validate_schema
from .types import SchemaValidationError, StepError, StepIterable


def run_steps(
    df: pd.DataFrame,
    steps: StepIterable,
    *,
    schema: DataFrameSchema | None = None,
    pipeline_version: str | None = None,
    logger=None,
) -> tuple[pd.DataFrame, logging_utils.PipelineRunMetrics]:
    """Execute ``steps`` sequentially and return the DataFrame and metrics."""

    log = logger or logging_utils.get_logger()
    current = ensure_dataframe(df, copy=True)
    if pipeline_version is not None:
        current.attrs["pipeline_version"] = pipeline_version

    started_at = datetime.now(timezone.utc).isoformat()
    pipeline_timer = perf_counter()
    metrics = logging_utils.PipelineRunMetrics(
        pipeline_version=pipeline_version,
        started_at=started_at,
        input_rows=current.shape[0],
        input_columns=current.shape[1],
    )

    for step in steps:
        try:
            next_df, step_metrics = logging_utils.execute_step(
                step.name,
                step.func,
                clone_dataframe(current),
                logger=log,
            )
        except SchemaValidationError:
            raise
        except StepError:
            raise
        except TypeError as exc:  # pragma: no cover - defensive
            raise StepError(step.name, str(exc), cause=exc) from exc
        except Exception as exc:  # pragma: no cover - defensive
            raise StepError(step.name, str(exc), cause=exc) from exc

        next_df = ensure_dataframe(next_df, copy=True)
        if pipeline_version is not None:
            next_df.attrs["pipeline_version"] = pipeline_version
        metrics.steps.append(step_metrics)
        current = next_df

    if schema is not None:
        schema_name = schema.__class__.__name__
        validation_start = perf_counter()
        validation_started_at = datetime.now(timezone.utc).isoformat()
        log.info("Validating final schema (%s)", schema_name)
        current = validate_schema(current, schema, context="pipeline")
        validation_duration = perf_counter() - validation_start
        validation_completed_at = datetime.now(timezone.utc).isoformat()
        if pipeline_version is not None:
            current.attrs["pipeline_version"] = pipeline_version
        metrics.validation = logging_utils.ValidationMetrics(
            schema=schema_name,
            started_at=validation_started_at,
            completed_at=validation_completed_at,
            duration_s=validation_duration,
        )
        log.info(
            "Schema validation completed | schema=%s | duration=%.3fs | rows=%s | cols=%s",
            schema_name,
            validation_duration,
            current.shape[0],
            current.shape[1],
        )

    metrics.finalize(
        output_rows=current.shape[0],
        output_columns=current.shape[1],
        duration_s=perf_counter() - pipeline_timer,
    )
    return current, metrics


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
            output=str(output_path),
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
            output=str(output_path),
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
            output=str(output_path),
            error=str(exc),
        )
        return None, None

    report_dir = Path(output_dir)
    report_path = report_dir / f"{table}.postprocess.report.json"
    payload = logging_utils.build_report_payload(
        table=table,
        metrics=metrics,
        output_path=str(output_path),
        extras=report_extras,
    )
    logging_utils.dump_report(report_path, payload)
    return metrics, report_path


__all__ = ["collect_postprocess_metrics", "infer_pipeline_version", "run_steps"]
