"""Execution helpers for postprocessing pipelines."""
from __future__ import annotations

import time
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd

from . import logging as logging_utils
from .config import ConfiguredStep
from .io import clone_dataframe, ensure_dataframe
from .schema import DataFrameSchema, validate_schema
from .types import SchemaValidationError, StepDefinition, StepError, StepFn

__all__ = [
    "PipelineRunMetadata",
    "PipelineStepStats",
    "SchemaValidationStats",
    "run_steps",
]


@dataclass(slots=True)
class SchemaValidationStats:
    """Metadata describing a schema validation pass."""

    name: str
    duration_seconds: float
    error: str | None = None

    def as_dict(self) -> dict[str, Any]:
        """Serialize the validation statistics into a plain mapping."""

        return {
            "name": self.name,
            "duration_seconds": self.duration_seconds,
            "error": self.error,
        }


@dataclass(slots=True)
class PipelineStepStats:
    """Capture execution statistics for a single pipeline step."""

    name: str
    description: str | None
    params: Mapping[str, Any]
    duration_seconds: float
    rows_before: int
    rows_after: int
    rows_delta: int
    columns_before: int
    columns_after: int
    columns_delta: int
    added_columns: tuple[str, ...]
    removed_columns: tuple[str, ...]
    error: str | None = None

    def as_dict(self) -> dict[str, Any]:
        """Serialize the collected statistics into a mapping."""

        return {
            "name": self.name,
            "description": self.description,
            "params": dict(self.params),
            "duration_seconds": self.duration_seconds,
            "rows_before": self.rows_before,
            "rows_after": self.rows_after,
            "rows_delta": self.rows_delta,
            "columns_before": self.columns_before,
            "columns_after": self.columns_after,
            "columns_delta": self.columns_delta,
            "added_columns": list(self.added_columns),
            "removed_columns": list(self.removed_columns),
            "error": self.error,
        }


@dataclass(slots=True)
class PipelineRunMetadata:
    """Aggregate metadata captured while executing a pipeline."""

    pipeline_version: str | None
    total_duration_seconds: float
    input_rows: int
    input_columns: int
    output_rows: int
    output_columns: int
    steps: tuple[PipelineStepStats, ...]
    pre_validation: SchemaValidationStats | None = None
    post_validation: SchemaValidationStats | None = None

    def as_dict(self) -> dict[str, Any]:
        """Represent the metadata using JSON-serialisable structures."""

        return {
            "pipeline_version": self.pipeline_version,
            "total_duration_seconds": self.total_duration_seconds,
            "input_rows": self.input_rows,
            "input_columns": self.input_columns,
            "output_rows": self.output_rows,
            "output_columns": self.output_columns,
            "steps": [step.as_dict() for step in self.steps],
            "pre_validation": (
                None if self.pre_validation is None else self.pre_validation.as_dict()
            ),
            "post_validation": (
                None if self.post_validation is None else self.post_validation.as_dict()
            ),
        }


StepLike = (
    ConfiguredStep
    | StepDefinition
    | tuple[StepFn, Mapping[str, Any]]
    | StepFn
)


@dataclass(frozen=True)
class _StepDescriptor:
    """Normalised view of a pipeline step definition."""

    name: str
    func: StepFn
    params: Mapping[str, Any]
    description: str | None


def _coerce_step(step: StepLike) -> _StepDescriptor:
    """Return a :class:`_StepDescriptor` for ``step``."""

    if isinstance(step, ConfiguredStep):
        definition = step.definition
        description = definition.description or step.description
        return _StepDescriptor(
            name=definition.name,
            func=definition.func,
            params=dict(step.params),
            description=description,
        )

    if isinstance(step, StepDefinition):
        return _StepDescriptor(
            name=step.name,
            func=step.func,
            params={},
            description=step.description,
        )

    if hasattr(step, "definition") and hasattr(step, "params"):
        definition = getattr(step, "definition")
        params = getattr(step, "params")
        if isinstance(definition, StepDefinition) and callable(definition.func):
            description = definition.description
            return _StepDescriptor(
                name=definition.name,
                func=definition.func,
                params=dict(params or {}),
                description=description,
            )

    if hasattr(step, "name") and hasattr(step, "func"):
        func = getattr(step, "func")
        if callable(func):
            params = getattr(step, "params", {}) or {}
            description = getattr(step, "description", None)
            name = getattr(step, "name", None) or getattr(func, "__name__", repr(func))
            return _StepDescriptor(
                name=name,
                func=func,
                params=dict(params),
                description=description,
            )

    if isinstance(step, tuple):
        if len(step) != 2:
            raise TypeError("step tuples must be expressed as (callable, params)")
        func, params = step
        if not callable(func):
            raise TypeError("step tuple must start with a callable")
        if params is None:
            params = {}
        if not isinstance(params, Mapping):
            raise TypeError("step tuple parameters must be a mapping")
        name = getattr(func, "__name__", repr(func))
        return _StepDescriptor(
            name=name,
            func=func,
            params=dict(params),
            description=None,
        )

    if callable(step):
        name = getattr(step, "__name__", repr(step))
        return _StepDescriptor(name=name, func=step, params={}, description=None)

    raise TypeError(f"Unsupported step specification: {step!r}")


def _build_step_stats(
    descriptor: _StepDescriptor,
    *,
    duration: float,
    before_shape: Sequence[int],
    before_columns: Sequence[str],
    result: pd.DataFrame | None,
    error: str | None,
) -> PipelineStepStats:
    if result is None:
        rows_after, cols_after = before_shape
        added_columns: tuple[str, ...] = ()
        removed_columns: tuple[str, ...] = ()
    else:
        rows_after, cols_after = result.shape
        added_columns = tuple(col for col in result.columns if col not in before_columns)
        removed_columns = tuple(col for col in before_columns if col not in result.columns)

    return PipelineStepStats(
        name=descriptor.name,
        description=descriptor.description,
        params=MappingProxyType(dict(descriptor.params)),
        duration_seconds=duration,
        rows_before=before_shape[0],
        rows_after=rows_after,
        rows_delta=rows_after - before_shape[0],
        columns_before=before_shape[1],
        columns_after=cols_after,
        columns_delta=cols_after - before_shape[1],
        added_columns=added_columns,
        removed_columns=removed_columns,
        error=error,
    )


def _validate_schema(
    df: pd.DataFrame,
    schema: DataFrameSchema,
    *,
    context: str,
    logger,
) -> tuple[pd.DataFrame, SchemaValidationStats]:
    schema_name = schema.__class__.__name__
    logger.info("Validating %s schema (%s)", context, schema_name)
    start = time.perf_counter()
    try:
        validated = validate_schema(df, schema, context=context)
    except SchemaValidationError as exc:
        duration = time.perf_counter() - start
        stats = SchemaValidationStats(schema_name, duration, str(exc))
        logger.error(
            "%s schema validation failed after %.4fs: %s",
            context,
            duration,
            exc,
        )
        setattr(exc, "validation_stats", stats)
        raise
    else:
        duration = time.perf_counter() - start
        stats = SchemaValidationStats(schema_name, duration, None)
        logger.info(
            "%s schema validation succeeded in %.4fs",
            context,
            duration,
        )
        return ensure_dataframe(validated, copy=True), stats


def run_steps(
    df: pd.DataFrame,
    steps: Iterable[StepLike],
    *,
    pre_schema: DataFrameSchema | None = None,
    post_schema: DataFrameSchema | None = None,
    pipeline_version: str | None = None,
    logger=None,
) -> tuple[pd.DataFrame, PipelineRunMetadata]:
    """Execute ``steps`` sequentially returning the processed DataFrame and metadata."""

    log = logger or logging_utils.get_logger()
    current = ensure_dataframe(df, copy=True)
    current = clone_dataframe(current)

    if pipeline_version is not None:
        current.attrs["pipeline_version"] = pipeline_version

    initial_rows, initial_columns = current.shape
    total_start = time.perf_counter()
    step_stats: list[PipelineStepStats] = []

    pre_validation: SchemaValidationStats | None = None
    if pre_schema is not None:
        try:
            current, pre_validation = _validate_schema(
                current, pre_schema, context="pipeline_pre", logger=log
            )
        except SchemaValidationError:
            if pipeline_version is not None:
                current.attrs["pipeline_version"] = pipeline_version
            raise
        if pipeline_version is not None:
            current.attrs["pipeline_version"] = pipeline_version

    for raw_step in steps:
        descriptor = _coerce_step(raw_step)
        log.info("Starting step %s", descriptor.name)
        before_columns = list(current.columns)
        before_shape = current.shape
        step_start = time.perf_counter()

        try:
            produced = descriptor.func(
                clone_dataframe(current), **dict(descriptor.params)
            )
        except SchemaValidationError as exc:
            duration = time.perf_counter() - step_start
            stats = _build_step_stats(
                descriptor,
                duration=duration,
                before_shape=before_shape,
                before_columns=before_columns,
                result=None,
                error=str(exc),
            )
            step_stats.append(stats)
            log.error("Step %s failed after %.4fs: %s", descriptor.name, duration, exc)
            raise
        except StepError as exc:
            duration = time.perf_counter() - step_start
            stats = _build_step_stats(
                descriptor,
                duration=duration,
                before_shape=before_shape,
                before_columns=before_columns,
                result=None,
                error=str(exc),
            )
            step_stats.append(stats)
            log.error("Step %s failed after %.4fs: %s", descriptor.name, duration, exc)
            raise
        except Exception as exc:  # pragma: no cover - defensive guard
            duration = time.perf_counter() - step_start
            message = f"{exc.__class__.__name__}: {exc}"
            stats = _build_step_stats(
                descriptor,
                duration=duration,
                before_shape=before_shape,
                before_columns=before_columns,
                result=None,
                error=message,
            )
            step_stats.append(stats)
            log.error("Step %s failed after %.4fs: %s", descriptor.name, duration, message)
            raise StepError(descriptor.name, message, cause=exc) from exc

        duration = time.perf_counter() - step_start

        if not isinstance(produced, pd.DataFrame):
            message = "step did not return a pandas.DataFrame"
            stats = _build_step_stats(
                descriptor,
                duration=duration,
                before_shape=before_shape,
                before_columns=before_columns,
                result=None,
                error=message,
            )
            step_stats.append(stats)
            log.error("Step %s failed after %.4fs: %s", descriptor.name, duration, message)
            raise StepError(descriptor.name, message)

        next_df = ensure_dataframe(produced, copy=True)
        if pipeline_version is not None:
            next_df.attrs["pipeline_version"] = pipeline_version

        stats = _build_step_stats(
            descriptor,
            duration=duration,
            before_shape=before_shape,
            before_columns=before_columns,
            result=next_df,
            error=None,
        )
        log.info(
            "Completed step %s in %.4fs | rows=%s->%s (%+d) | cols=%s->%s (%+d)",
            descriptor.name,
            duration,
            before_shape[0],
            next_df.shape[0],
            stats.rows_delta,
            before_shape[1],
            next_df.shape[1],
            stats.columns_delta,
        )
        step_stats.append(stats)
        current = next_df

    post_validation: SchemaValidationStats | None = None
    if post_schema is not None:
        try:
            current, post_validation = _validate_schema(
                current, post_schema, context="pipeline_post", logger=log
            )
        except SchemaValidationError:
            if pipeline_version is not None:
                current.attrs["pipeline_version"] = pipeline_version
            raise
        if pipeline_version is not None:
            current.attrs["pipeline_version"] = pipeline_version

    total_duration = time.perf_counter() - total_start
    output_rows, output_columns = current.shape
    effective_version = pipeline_version or current.attrs.get("pipeline_version")

    metadata = PipelineRunMetadata(
        pipeline_version=effective_version,
        total_duration_seconds=total_duration,
        input_rows=initial_rows,
        input_columns=initial_columns,
        output_rows=output_rows,
        output_columns=output_columns,
        steps=tuple(step_stats),
        pre_validation=pre_validation,
        post_validation=post_validation,
    )

    return current, metadata

