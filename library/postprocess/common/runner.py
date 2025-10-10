"""Orchestration helpers for executing post-processing pipelines."""

from __future__ import annotations

import inspect
import re
from collections.abc import Iterable, Mapping
from datetime import UTC, datetime
from time import perf_counter
from typing import Any

import pandas as pd

from . import logging as logging_utils
from .io import clone_dataframe, ensure_dataframe
from .schema import DataFrameSchema, validate_schema
from .types import SchemaValidationError, StepDefinition, StepError

StepTuple = tuple[Any, ...]


def _describe_dtypes(frame: pd.DataFrame) -> dict[str, str]:
    return {column: str(dtype) for column, dtype in frame.dtypes.items()}


def _now_iso() -> str:
    return datetime.now(UTC).isoformat()


def _coerce_step(entry: StepDefinition | StepTuple, index: int) -> StepDefinition:
    if isinstance(entry, StepDefinition):
        return entry

    if isinstance(entry, tuple):
        if len(entry) == 2:
            func, params = entry
            name = getattr(func, "__name__", f"step_{index}")
        elif len(entry) == 3:
            name, func, params = entry
        else:  # pragma: no cover - defensive guard
            raise TypeError(
                "Step tuples must be of the form (callable, params) or (name, callable, params)"
            )

        if not callable(func):
            raise TypeError(f"Step #{index} is not callable: {func!r}")
        if params is None:
            params = {}
        if not isinstance(params, Mapping):
            raise TypeError(
                f"Step #{index} parameters must be a mapping, got {type(params)!r}"
            )

        return StepDefinition(name=name, func=func, params=params)

    raise TypeError(f"Unsupported step entry at position {index}: {entry!r}")


def _normalise_steps(
    steps: Iterable[StepDefinition | StepTuple],
) -> tuple[StepDefinition, ...]:
    return tuple(_coerce_step(entry, index) for index, entry in enumerate(steps))


_UNEXPECTED_KEYWORD_PATTERN = re.compile(
    r"got an unexpected keyword argument ['\"](?P<name>[^'\"]+)['\"]"
)


def _extract_unexpected_keyword(error: TypeError) -> str | None:
    """Return the keyword name when ``error`` signals an unexpected argument."""

    message = str(error)
    match = _UNEXPECTED_KEYWORD_PATTERN.search(message)
    if match is not None:
        return match.group("name")
    return None


def _record_error(
    *,
    metrics: logging_utils.PipelineRunMetrics,
    step_name: str,
    started_at: str,
    start_clock: float,
    input_rows: int,
    input_columns: int,
    params: Mapping[str, Any],
    error: BaseException,
) -> None:
    metrics.steps.append(
        logging_utils.StepMetrics(
            name=step_name,
            started_at=started_at,
            completed_at=_now_iso(),
            duration_s=perf_counter() - start_clock,
            input_rows=input_rows,
            input_columns=input_columns,
            output_rows=input_rows,
            output_columns=input_columns,
            schema_diff=logging_utils.SchemaDiff(),
            parameters=dict(params),
            error=logging_utils.StepErrorInfo(
                type=error.__class__.__name__,
                message=str(error),
            ),
        )
    )


def _prepare_step_arguments(
    func: Any,
    params: Mapping[str, Any],
    *,
    logger,
    step_name: str,
) -> dict[str, Any]:
    """Return a mapping limited to arguments accepted by ``func``."""

    if not params:
        return {}

    try:
        signature = inspect.signature(func)
    except (TypeError, ValueError):  # pragma: no cover - fallback for builtins
        return dict(params)

    has_var_kwargs = any(
        parameter.kind is inspect.Parameter.VAR_KEYWORD
        for parameter in signature.parameters.values()
    )

    accepted: dict[str, Any] = {}
    ignored: list[str] = []

    for key, value in params.items():
        if key in signature.parameters or has_var_kwargs:
            accepted[key] = value
        else:
            ignored.append(key)

    if ignored:
        logger.warning(
            "Step %s ignoring unsupported parameters: %s",
            step_name,
            ", ".join(sorted(ignored)),
        )

    return accepted


def _identify_unsupported_kwargs(
    func: Any, params: Mapping[str, Any]
) -> tuple[str, ...]:
    """Return a tuple of keyword names not accepted by ``func``."""

    if not params:
        return ()

    try:
        signature = inspect.signature(func)
    except (TypeError, ValueError):  # pragma: no cover - defensive
        return ()

    if any(
        parameter.kind is inspect.Parameter.VAR_KEYWORD
        for parameter in signature.parameters.values()
    ):
        return ()

    unsupported: list[str] = []
    for key in params:
        if key not in signature.parameters:
            unsupported.append(key)
    return tuple(unsupported)


def run_steps(
    df: pd.DataFrame,
    steps: Iterable[StepDefinition | StepTuple],
    *,
    pipeline_version: str | None = None,
    pre_schema: DataFrameSchema | None = None,
    post_schema: DataFrameSchema | None = None,
    logger=None,
) -> tuple[pd.DataFrame, logging_utils.PipelineRunMetrics]:
    """Execute ``steps`` sequentially and return the transformed frame and metadata."""

    log = logger or logging_utils.get_logger()
    current = ensure_dataframe(df, copy=True)
    if pipeline_version is not None:
        current.attrs["pipeline_version"] = pipeline_version

    normalised_steps = _normalise_steps(steps)

    pipeline_started_at = _now_iso()
    pipeline_timer = perf_counter()
    metrics = logging_utils.PipelineRunMetrics(
        pipeline_version=pipeline_version,
        started_at=pipeline_started_at,
        input_rows=current.shape[0],
        input_columns=current.shape[1],
    )

    if pre_schema is not None:
        schema_name = pre_schema.__class__.__name__
        log.info("Validating input schema (%s)", schema_name)
        current = validate_schema(current, pre_schema, context="pre_pipeline")
        if pipeline_version is not None:
            current.attrs["pipeline_version"] = pipeline_version

    for index, step in enumerate(normalised_steps):
        params = dict(step.params)
        call_params = _prepare_step_arguments(
            step.func,
            params,
            logger=log,
            step_name=step.name,
        )
        call_params = dict(call_params)
        frame = clone_dataframe(current)
        step_started_at = _now_iso()
        step_clock = perf_counter()
        input_rows, input_columns = frame.shape
        before_columns = list(frame.columns)
        before_dtypes = _describe_dtypes(frame)

        log.info("Starting step %s", step.name)

        try:
            while True:
                try:
                    result = step.func(frame, **call_params)
                    break
                except TypeError as exc:
                    unexpected_keyword = _extract_unexpected_keyword(exc)
                    if (
                        unexpected_keyword is None
                        or unexpected_keyword not in call_params
                    ):
                        fallback = _identify_unsupported_kwargs(step.func, call_params)
                        if not fallback:
                            raise
                        log.warning(
                            "Step %s retrying without unsupported parameters: %s",
                            step.name,
                            ", ".join(sorted(fallback)),
                        )
                        for key in fallback:
                            call_params.pop(key, None)
                        continue

                    log.warning(
                        "Step %s retrying without unsupported parameter: %s",
                        step.name,
                        unexpected_keyword,
                    )
                    call_params.pop(unexpected_keyword, None)

            # normal execution path continues below once the loop breaks
        except SchemaValidationError as exc:
            log.exception("Schema validation failed during step %s", step.name)
            _record_error(
                metrics=metrics,
                step_name=step.name,
                started_at=step_started_at,
                start_clock=step_clock,
                input_rows=input_rows,
                input_columns=input_columns,
                params=params,
                error=exc,
            )
            raise
        except StepError as exc:
            log.exception("Step %s raised StepError", step.name)
            _record_error(
                metrics=metrics,
                step_name=step.name,
                started_at=step_started_at,
                start_clock=step_clock,
                input_rows=input_rows,
                input_columns=input_columns,
                params=params,
                error=exc,
            )
            raise
        except TypeError as exc:  # pragma: no cover - defensive guard
            log.exception("Step %s encountered a TypeError", step.name)
            _record_error(
                metrics=metrics,
                step_name=step.name,
                started_at=step_started_at,
                start_clock=step_clock,
                input_rows=input_rows,
                input_columns=input_columns,
                params=params,
                error=exc,
            )
            raise StepError(step.name, str(exc), cause=exc) from exc
        except Exception as exc:  # pragma: no cover - defensive guard
            log.exception("Step %s failed with an unexpected error", step.name)
            _record_error(
                metrics=metrics,
                step_name=step.name,
                started_at=step_started_at,
                start_clock=step_clock,
                input_rows=input_rows,
                input_columns=input_columns,
                params=params,
                error=exc,
            )
            raise StepError(step.name, str(exc), cause=exc) from exc

        next_df = ensure_dataframe(result, copy=True)
        completed_at = _now_iso()
        duration = perf_counter() - step_clock
        output_rows, output_columns = next_df.shape
        after_columns = list(next_df.columns)
        after_dtypes = _describe_dtypes(next_df)
        diff = logging_utils.compute_schema_diff(
            before_columns,
            before_dtypes,
            after_columns,
            after_dtypes,
        )

        log.info(
            "Completed step %s | duration=%.3fs | rows=%s->%s | cols=%s->%s | added=%s | removed=%s | types=%s",
            step.name,
            duration,
            input_rows,
            output_rows,
            input_columns,
            output_columns,
            ",".join(diff.added) or "-",
            ",".join(diff.removed) or "-",
            ", ".join(
                f"{column}: {change['from']} -> {change['to']}"
                for column, change in diff.type_changes.items()
            )
            or "-",
        )

        if pipeline_version is not None:
            next_df.attrs["pipeline_version"] = pipeline_version

        metrics.steps.append(
            logging_utils.StepMetrics(
                name=step.name,
                started_at=step_started_at,
                completed_at=completed_at,
                duration_s=duration,
                input_rows=input_rows,
                input_columns=input_columns,
                output_rows=output_rows,
                output_columns=output_columns,
                schema_diff=diff,
                parameters=params,
            )
        )

        current = next_df

    if post_schema is not None:
        schema_name = post_schema.__class__.__name__
        validation_started_at = _now_iso()
        validation_clock = perf_counter()
        log.info("Validating final schema (%s)", schema_name)
        current = validate_schema(current, post_schema, context="pipeline")
        validation_duration = perf_counter() - validation_clock
        validation_completed_at = _now_iso()
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

    if pipeline_version is not None:
        current.attrs["pipeline_version"] = pipeline_version

    return current, metrics


__all__ = ["run_steps"]
