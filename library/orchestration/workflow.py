"""Helpers to execute prepared pipeline steps without invoking CLI parsers."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from library.pipelines.registry import PipelineStep

if TYPE_CHECKING:  # pragma: no cover - imports used only for typing
    from scripts.get_data import PipelineRunConfig


@dataclass(frozen=True)
class StepExecutionResult:
    """Summarize the outcome of invoking a pipeline step."""

    exit_code: int
    executed: bool
    status: str
    reason: str | None = None


StepCallable = Callable[["PipelineRunConfig", Path, Path, Path], StepExecutionResult]


@dataclass(frozen=True)
class PreparedPipelineStep:
    """Couple :class:`PipelineStep` metadata with an execution callable."""

    step: PipelineStep
    invoke: StepCallable


@dataclass(frozen=True)
class WorkflowStepResult:
    """Container describing the outcome of a prepared pipeline step."""

    step: PipelineStep
    input_path: Path
    final_output: Path
    working_output: Path
    result: StepExecutionResult


def temporary_output_path(output_path: Path) -> Path:
    """Return the path used for intermediate artefacts of ``output_path``."""

    return output_path.with_name(f".{output_path.name}.tmp")


def execute_workflow(
    cfg: PipelineRunConfig,
    steps: Sequence[PreparedPipelineStep],
) -> Iterator[WorkflowStepResult]:
    """Execute ``steps`` sequentially and yield their results.

    The provided callables are expected to handle the step invocation without
    triggering command-line parsing. ``execute_workflow`` computes the
    filesystem paths that each callable would typically derive from CLI
    arguments and passes them directly, ensuring deterministic execution order
    and consistent input/output handling.
    """

    for prepared in steps:
        step = prepared.step
        input_path = step.required_input(cfg)
        final_output = step.expected_output(cfg)
        working_output = temporary_output_path(final_output)
        result = prepared.invoke(cfg, input_path, final_output, working_output)
        yield WorkflowStepResult(
            step=step,
            input_path=input_path,
            final_output=final_output,
            working_output=working_output,
            result=result,
        )
        if result.exit_code != 0:
            break


__all__ = [
    "PreparedPipelineStep",
    "StepCallable",
    "StepExecutionResult",
    "WorkflowStepResult",
    "execute_workflow",
    "temporary_output_path",
]
