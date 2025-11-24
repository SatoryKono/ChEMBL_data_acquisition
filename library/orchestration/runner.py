"""Helpers for preparing and executing pipeline steps.

This module centralises common logic for translating :class:`PipelineStep`
definitions into concrete filesystem paths and for invoking prepared steps in
a deterministic sequence. Callers can reuse these helpers to avoid duplicating
path preparation or artefact resolution when orchestrating pipelines outside of
the CLI entrypoints.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence, TYPE_CHECKING

from library.orchestration.workflow import (
    PreparedPipelineStep,
    WorkflowStepResult,
    temporary_output_path,
)
from library.pipelines.registry import PipelineStep

if TYPE_CHECKING:  # pragma: no cover - imports only needed for typing
    from library.cli.commands.get_data import PipelineRunConfig

__all__ = [
    "StepPaths",
    "failure_sentinel_path",
    "prepare_step_paths",
    "resolve_consumed_artifact_path",
    "run_prepared_steps",
]


@dataclass(frozen=True)
class StepPaths:
    """Filesystem locations associated with a pipeline step."""

    input_path: Path
    final_output: Path
    working_output: Path
    sentinel_path: Path


def failure_sentinel_path(output_path: Path) -> Path:
    """Return the sentinel path recorded when a pipeline step fails."""

    return output_path.with_name(f"{output_path.name}.failed")


def resolve_consumed_artifact_path(cfg: "PipelineRunConfig", artefact: str) -> Path:
    """Return the filesystem path associated with a consumed artefact name."""

    candidate = Path(artefact)
    if candidate.is_absolute():
        return candidate
    if candidate.suffix:
        return cfg.output_dir / candidate
    return cfg.output_dir / f"output.{candidate.name}_{cfg.date_prefix}.csv"


def prepare_step_paths(step: PipelineStep, cfg: "PipelineRunConfig") -> StepPaths:
    """Compute input/output paths for ``step`` using ``cfg``."""

    final_output = step.expected_output(cfg)
    working_output = temporary_output_path(final_output)
    return StepPaths(
        input_path=step.required_input(cfg),
        final_output=final_output,
        working_output=working_output,
        sentinel_path=failure_sentinel_path(final_output),
    )


def run_prepared_steps(
    cfg: "PipelineRunConfig", steps: Sequence[PreparedPipelineStep]
) -> Iterator[WorkflowStepResult]:
    """Execute ``steps`` sequentially using prepared path arguments."""

    for prepared in steps:
        paths = prepare_step_paths(prepared.step, cfg)
        result = prepared.invoke(
            cfg, paths.input_path, paths.final_output, paths.working_output
        )
        yield WorkflowStepResult(
            step=prepared.step,
            input_path=paths.input_path,
            final_output=paths.final_output,
            working_output=paths.working_output,
            result=result,
        )
        if result.exit_code != 0:
            break
