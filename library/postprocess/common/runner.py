"""Execution utilities for postprocessing pipelines."""
from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from types import MappingProxyType
from typing import Any, Literal, Mapping

import pandas as pd

from . import logging as logging_utils
from .io import clone_dataframe, ensure_dataframe
from .schema import DataFrameSchema, validate_schema
from .types import (
    SchemaValidationError,
    StepDefinition,
    StepError,
    StepIterable,
    StepLike,
)


@dataclass(frozen=True)
class SchemaCheckReport:
    """Summary of a schema validation stage."""

    schema: str
    duration_seconds: float
    status: Literal["success", "error"]
    error: str | None = None


@dataclass(frozen=True)
class StepReport:
    """Execution statistics for a single pipeline step."""

    name: str
    parameters: Mapping[str, Any]
    duration_seconds: float
    rows_before: int
    rows_after: int | None
    row_delta: int | None
    columns_before: tuple[str, ...]
    columns_after: tuple[str, ...] | None
    column_delta: int | None
    added_columns: tuple[str, ...]
    removed_columns: tuple[str, ...]
    error: str | None = None


@dataclass(frozen=True)
class RunnerMetadata:
    """Aggregated metadata about a pipeline run."""

    pipeline_version: str | None
    steps: tuple[StepReport, ...]
    pre_schema: SchemaCheckReport | None
    post_schema: SchemaCheckReport | None


@dataclass(frozen=True)
class RunnerResult:
    """Return value of :func:`run_steps`."""

    frame: pd.DataFrame
    metadata: RunnerMetadata


def _normalise_step(step: StepLike) -> StepDefinition:
    if isinstance(step, StepDefinition):
        return step

    if callable(step):
        name = getattr(step, "__name__", "step")
        return StepDefinition(name, step)

    if isinstance(step, tuple):
        if len(step) == 1:
            (func,) = step
            if not callable(func):
                raise TypeError("Step tuple must contain callables")
            name = getattr(func, "__name__", "step")
            return StepDefinition(name, func)
        if len(step) == 2:
            first, second = step
            if isinstance(first, str) and callable(second):
                return StepDefinition(first, second)
            if callable(first):
                func = first
                params = second
                if params is not None and not isinstance(params, Mapping):
                    raise TypeError("Step parameters must be provided as a mapping")
                name = getattr(func, "__name__", "step")
                return StepDefinition(name, func, parameters=params)
        if len(step) == 3:
            name, func, params = step
            if not isinstance(name, str) or not callable(func):
                raise TypeError("Step tuple must be (name, callable, params)")
            if params is not None and not isinstance(params, Mapping):
                raise TypeError("Step parameters must be provided as a mapping")
            return StepDefinition(name, func, parameters=params)

    raise TypeError(f"Unsupported step specification: {step!r}")


def _coerce_steps(steps: StepIterable) -> tuple[StepDefinition, ...]:
    return tuple(_normalise_step(step) for step in steps)


def _sorted_parameters(params: Mapping[str, Any]) -> Mapping[str, Any]:
    if not params:
        return MappingProxyType({})
    ordered = {key: params[key] for key in sorted(params)}
    return MappingProxyType(ordered)


def run_steps(
    df: pd.DataFrame,
    steps: StepIterable,
    *,
    pre_schema: DataFrameSchema | None = None,
    post_schema: DataFrameSchema | None = None,
    pipeline_version: str | None = None,
    logger=None,
) -> RunnerResult:
    """Execute ``steps`` sequentially returning the processed DataFrame and metadata."""

    log = logger or logging_utils.get_logger()
    current = ensure_dataframe(df, copy=True)
    if pipeline_version is not None:
        current.attrs["pipeline_version"] = pipeline_version

    step_reports: list[StepReport] = []
    pre_schema_report: SchemaCheckReport | None = None
    post_schema_report: SchemaCheckReport | None = None

    if pre_schema is not None:
        log.info("Validating input schema (%s)", pre_schema.__class__.__name__)
        start = perf_counter()
        try:
            current = validate_schema(current, pre_schema, context="pipeline.pre_schema")
        except SchemaValidationError as exc:
            duration = perf_counter() - start
            pre_schema_report = SchemaCheckReport(
                schema=pre_schema.__class__.__name__,
                duration_seconds=duration,
                status="error",
                error=f"{exc.__class__.__name__}: {exc}",
            )
            raise
        else:
            duration = perf_counter() - start
            pre_schema_report = SchemaCheckReport(
                schema=pre_schema.__class__.__name__,
                duration_seconds=duration,
                status="success",
            )
            if pipeline_version is not None:
                current.attrs["pipeline_version"] = pipeline_version

    for step in _coerce_steps(steps):
        log.info("Starting step %s", step.name)
        before_columns: tuple[str, ...] = tuple(current.columns)
        before_rows = current.shape[0]
        start = perf_counter()
        params = dict(step.parameters)
        error_message: str | None = None
        try:
            next_df = step.func(clone_dataframe(current), **params)
        except SchemaValidationError as exc:
            error_message = f"{exc.__class__.__name__}: {exc}"
            duration = perf_counter() - start
            step_reports.append(
                StepReport(
                    name=step.name,
                    parameters=_sorted_parameters(params),
                    duration_seconds=duration,
                    rows_before=before_rows,
                    rows_after=None,
                    row_delta=None,
                    columns_before=before_columns,
                    columns_after=None,
                    column_delta=None,
                    added_columns=(),
                    removed_columns=(),
                    error=error_message,
                )
            )
            raise
        except StepError as exc:
            error_message = f"{exc.__class__.__name__}: {exc}"
            duration = perf_counter() - start
            step_reports.append(
                StepReport(
                    name=step.name,
                    parameters=_sorted_parameters(params),
                    duration_seconds=duration,
                    rows_before=before_rows,
                    rows_after=None,
                    row_delta=None,
                    columns_before=before_columns,
                    columns_after=None,
                    column_delta=None,
                    added_columns=(),
                    removed_columns=(),
                    error=error_message,
                )
            )
            raise
        except Exception as exc:  # pragma: no cover - defensive
            error_message = f"{exc.__class__.__name__}: {exc}"
            duration = perf_counter() - start
            step_reports.append(
                StepReport(
                    name=step.name,
                    parameters=_sorted_parameters(params),
                    duration_seconds=duration,
                    rows_before=before_rows,
                    rows_after=None,
                    row_delta=None,
                    columns_before=before_columns,
                    columns_after=None,
                    column_delta=None,
                    added_columns=(),
                    removed_columns=(),
                    error=error_message,
                )
            )
            raise StepError(step.name, str(exc), cause=exc) from exc

        duration = perf_counter() - start
        if not isinstance(next_df, pd.DataFrame):
            error_message = "Step did not return a pandas.DataFrame"
            step_reports.append(
                StepReport(
                    name=step.name,
                    parameters=_sorted_parameters(params),
                    duration_seconds=duration,
                    rows_before=before_rows,
                    rows_after=None,
                    row_delta=None,
                    columns_before=before_columns,
                    columns_after=None,
                    column_delta=None,
                    added_columns=(),
                    removed_columns=(),
                    error=error_message,
                )
            )
            raise StepError(step.name, error_message)

        next_df = ensure_dataframe(next_df, copy=True)
        if pipeline_version is not None:
            next_df.attrs["pipeline_version"] = pipeline_version

        after_columns: tuple[str, ...] = tuple(next_df.columns)
        after_rows = next_df.shape[0]
        added_columns = tuple(col for col in after_columns if col not in before_columns)
        removed_columns = tuple(col for col in before_columns if col not in after_columns)
        row_delta = after_rows - before_rows
        column_delta = len(after_columns) - len(before_columns)

        log.info(
            "Completed step %s | rows=%s->%s | cols=%s->%s | +%s | -%s",
            step.name,
            before_rows,
            after_rows,
            len(before_columns),
            len(after_columns),
            ",".join(added_columns) or "-",
            ",".join(removed_columns) or "-",
        )

        step_reports.append(
            StepReport(
                name=step.name,
                parameters=_sorted_parameters(params),
                duration_seconds=duration,
                rows_before=before_rows,
                rows_after=after_rows,
                row_delta=row_delta,
                columns_before=before_columns,
                columns_after=after_columns,
                column_delta=column_delta,
                added_columns=added_columns,
                removed_columns=removed_columns,
                error=None,
            )
        )

        current = next_df

    if post_schema is not None:
        log.info("Validating output schema (%s)", post_schema.__class__.__name__)
        start = perf_counter()
        try:
            current = validate_schema(current, post_schema, context="pipeline.post_schema")
        except SchemaValidationError as exc:
            duration = perf_counter() - start
            post_schema_report = SchemaCheckReport(
                schema=post_schema.__class__.__name__,
                duration_seconds=duration,
                status="error",
                error=f"{exc.__class__.__name__}: {exc}",
            )
            raise
        else:
            duration = perf_counter() - start
            post_schema_report = SchemaCheckReport(
                schema=post_schema.__class__.__name__,
                duration_seconds=duration,
                status="success",
            )
            if pipeline_version is not None:
                current.attrs["pipeline_version"] = pipeline_version

    metadata = RunnerMetadata(
        pipeline_version=pipeline_version,
        steps=tuple(step_reports),
        pre_schema=pre_schema_report,
        post_schema=post_schema_report,
    )

    return RunnerResult(frame=current, metadata=metadata)


__all__ = [
    "RunnerMetadata",
    "RunnerResult",
    "SchemaCheckReport",
    "StepReport",
    "run_steps",
]
