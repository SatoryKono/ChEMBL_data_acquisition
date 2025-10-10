"""Logging utilities for postprocessing pipelines."""

from __future__ import annotations

import json
import logging
import os
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from time import perf_counter
from typing import Any

import pandas as pd

_DEFAULT_LOGGER_NAME = "chembl.postprocess"


def _now_iso() -> str:
    """Return an ISO formatted UTC timestamp."""

    return datetime.now(UTC).isoformat()


def _duration_seconds(start: float, end: float | None = None) -> float:
    """Return the duration in seconds between ``start`` and ``end``."""

    stop = perf_counter() if end is None else end
    return max(0.0, stop - start)


def _describe_dtypes(frame: pd.DataFrame) -> dict[str, str]:
    """Return a mapping of column name to dtype string for ``frame``."""

    return {column: str(dtype) for column, dtype in frame.dtypes.items()}


@dataclass(slots=True)
class SchemaDiff:
    """Representation of column-level schema mutations within a step."""

    added: list[str] = field(default_factory=list)
    removed: list[str] = field(default_factory=list)
    type_changes: dict[str, dict[str, str]] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serialisable representation of the diff."""

        return {
            "added": list(self.added),
            "removed": list(self.removed),
            "type_changes": {
                column: {"from": change["from"], "to": change["to"]}
                for column, change in self.type_changes.items()
            },
        }


@dataclass(slots=True)
class StepErrorInfo:
    """Structured representation of a step failure."""

    type: str
    message: str

    def to_dict(self) -> dict[str, Any]:
        return {"type": self.type, "message": self.message}


@dataclass(slots=True)
class StepMetrics:
    """Metrics captured for a single transformation step."""

    name: str
    started_at: str
    completed_at: str
    duration_s: float
    input_rows: int
    input_columns: int
    output_rows: int
    output_columns: int
    schema_diff: SchemaDiff = field(default_factory=SchemaDiff)
    parameters: dict[str, Any] = field(default_factory=dict)
    error: StepErrorInfo | None = None

    def to_dict(self) -> dict[str, Any]:
        """Return a dictionary representation suitable for JSON output."""

        return {
            "name": self.name,
            "started_at": self.started_at,
            "completed_at": self.completed_at,
            "duration_s": self.duration_s,
            "input": {"rows": self.input_rows, "columns": self.input_columns},
            "output": {"rows": self.output_rows, "columns": self.output_columns},
            "schema_diff": self.schema_diff.to_dict(),
            "parameters": dict(self.parameters),
            "error": self.error.to_dict() if self.error else None,
        }


@dataclass(slots=True)
class ValidationMetrics:
    """Metrics describing the terminal schema validation step."""

    schema: str
    started_at: str
    completed_at: str
    duration_s: float
    status: str = "passed"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "status": self.status,
            "started_at": self.started_at,
            "completed_at": self.completed_at,
            "duration_s": self.duration_s,
        }


@dataclass(slots=True)
class PipelineRunMetrics:
    """Aggregate metrics for the full postprocessing pipeline."""

    pipeline_version: str | None
    started_at: str
    input_rows: int
    input_columns: int
    steps: list[StepMetrics] = field(default_factory=list)
    validation: ValidationMetrics | None = None
    completed_at: str | None = None
    duration_s: float | None = None
    output_rows: int | None = None
    output_columns: int | None = None

    def finalize(
        self, *, output_rows: int, output_columns: int, duration_s: float
    ) -> None:
        """Record terminal metrics once the pipeline has finished."""

        self.output_rows = output_rows
        self.output_columns = output_columns
        self.completed_at = _now_iso()
        self.duration_s = duration_s

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serialisable view of the captured metrics."""

        return {
            "pipeline_version": self.pipeline_version,
            "started_at": self.started_at,
            "completed_at": self.completed_at,
            "duration_s": self.duration_s,
            "input": {"rows": self.input_rows, "columns": self.input_columns},
            "output": {
                "rows": self.output_rows,
                "columns": self.output_columns,
            },
            "steps": [step.to_dict() for step in self.steps],
            "validation": self.validation.to_dict() if self.validation else None,
        }

    def summary(self) -> dict[str, Any]:
        """Return a lightweight summary useful for structured logging."""

        return {
            "pipeline_version": self.pipeline_version,
            "steps": len(self.steps),
            "duration_s": self.duration_s,
            "rows": self.output_rows,
            "columns": self.output_columns,
        }


def get_logger(name: str | None = None) -> logging.Logger:
    """Return a configured logger for postprocessing steps.

    The logger defaults to :mod:`chembl.postprocess` and ensures that
    double-configuration is avoided when used inside notebooks or CLI tools.
    """

    logger_name = name or _DEFAULT_LOGGER_NAME
    logger = logging.getLogger(logger_name)

    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
            datefmt="%Y-%m-%dT%H:%M:%S",
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)

    logger.propagate = False
    return logger


def compute_schema_diff(
    before_columns: list[str],
    before_dtypes: Mapping[str, str],
    after_columns: list[str],
    after_dtypes: Mapping[str, str],
) -> SchemaDiff:
    """Return a :class:`SchemaDiff` describing changes between two frames."""

    added = [column for column in after_columns if column not in before_columns]
    removed = [column for column in before_columns if column not in after_columns]
    type_changes: dict[str, dict[str, str]] = {}
    for column in before_dtypes.keys() & after_dtypes.keys():
        before_type = before_dtypes[column]
        after_type = after_dtypes[column]
        if before_type != after_type:
            type_changes[column] = {"from": before_type, "to": after_type}
    return SchemaDiff(added=added, removed=removed, type_changes=type_changes)


def execute_step(
    step_name: str,
    func: Callable[[pd.DataFrame], pd.DataFrame],
    frame: pd.DataFrame,
    *,
    logger: logging.Logger,
    params: Mapping[str, Any] | None = None,
) -> tuple[pd.DataFrame, StepMetrics]:
    """Execute ``func`` capturing timing and schema mutations."""

    start_clock = perf_counter()
    started_at = _now_iso()
    before_columns = list(frame.columns)
    before_dtypes = _describe_dtypes(frame)
    input_rows, input_columns = frame.shape

    parameters = dict(params or {})
    logger.info("Starting step %s", step_name)
    result = func(frame, **parameters)

    completed_at = _now_iso()
    duration_s = _duration_seconds(start_clock)
    if not isinstance(result, pd.DataFrame):
        raise TypeError(f"Step '{step_name}' did not return a pandas.DataFrame")

    output_rows, output_columns = result.shape
    after_columns = list(result.columns)
    after_dtypes = _describe_dtypes(result)
    diff = compute_schema_diff(
        before_columns, before_dtypes, after_columns, after_dtypes
    )

    type_change_summary = ", ".join(
        f"{column}: {change['from']} -> {change['to']}"
        for column, change in diff.type_changes.items()
    )
    logger.info(
        "Completed step %s | duration=%.3fs | rows=%s->%s | cols=%s->%s | added=%s | removed=%s | types=%s",
        step_name,
        duration_s,
        input_rows,
        output_rows,
        input_columns,
        output_columns,
        ",".join(diff.added) or "-",
        ",".join(diff.removed) or "-",
        type_change_summary or "-",
    )

    metrics = StepMetrics(
        name=step_name,
        started_at=started_at,
        completed_at=completed_at,
        duration_s=duration_s,
        input_rows=input_rows,
        input_columns=input_columns,
        output_rows=output_rows,
        output_columns=output_columns,
        schema_diff=diff,
        parameters=parameters,
    )
    return result, metrics


def build_report_payload(
    table: str,
    metrics: PipelineRunMetrics,
    *,
    output_path: str | None = None,
    output_postprocessed: str | None = None,
    extras: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a JSON payload representing the collected metrics."""

    payload = {
        "table": table,
        "metrics": metrics.to_dict(),
    }

    # `output_path` is a legacy key consumed by existing tooling.
    # `output_postprocessed` is the new preferred field.
    # For backwards compatibility, always provide both fields if possible:
    effective_postprocessed = (
        output_postprocessed if output_postprocessed is not None else output_path
    )
    if effective_postprocessed is not None:
        payload["output_postprocessed"] = effective_postprocessed
    if effective_postprocessed is not None:
        payload["output_path"] = effective_postprocessed
    if extras:
        payload["extras"] = dict(extras)
    return payload


def dump_report(
    path: str | os.PathLike[str],
    payload: Mapping[str, Any],
    *,
    indent: int = 2,
) -> None:
    """Write ``payload`` as JSON to ``path``."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(payload, indent=indent, sort_keys=True), encoding="utf-8"
    )


__all__ = [
    "PipelineRunMetrics",
    "SchemaDiff",
    "StepErrorInfo",
    "StepMetrics",
    "ValidationMetrics",
    "build_report_payload",
    "compute_schema_diff",
    "dump_report",
    "execute_step",
    "get_logger",
]
