from __future__ import annotations

import argparse
import importlib
from collections.abc import Callable, Sequence
from types import ModuleType
from typing import cast

from library.common.log import logger
<<<<<<< HEAD
from library.common.logging_setup import Logger
=======
from library.config import Config
>>>>>>> origin/codex/fix-styling-baseline-in-ci-7cexye
from library.pipelines.activity.runner import (
    MIN_ACTIVITY_TIMEOUT,
    ActivityCommandOptions,
    resolve_activity_pipeline_hooks,
)
from library.pipelines.activity.runner import (
    run_activity_pipeline as _run_activity_pipeline,
)

from . import _run

<<<<<<< HEAD
def _sync_pipeline_logger(current_logger: Logger) -> None:
    """Align the shared pipeline logger with ``current_logger`` when possible."""
=======

def _sync_pipeline_logger(logger: object) -> None:
    """Align the shared pipeline logger with ``logger`` when possible."""
>>>>>>> origin/codex/fix-styling-baseline-in-ci-7cexye

    modules: Sequence[tuple[str, str]] = (
        ("library.common.log", "logger"),
        ("library.pipelines.activity.runner", "logger"),
    )
    for module_name, attribute in modules:
        try:
            module = importlib.import_module(module_name)
        except Exception:  # pragma: no cover - defensive
            continue
        module_typed = cast(ModuleType, module)
        if getattr(module_typed, attribute, None) is current_logger:
            continue
        try:
            setattr(module_typed, attribute, current_logger)
        except Exception:  # pragma: no cover - defensive
            continue


def run_activity_pipeline(
    cfg: Config,
    options: ActivityCommandOptions,
    *,
    runner: Callable[[Config, argparse.Namespace], int] | None = None,
    emit_completion_message: Callable[..., None] | None = None,
) -> int:
    """Execute the activity pipeline with orchestrator semantics."""

    if runner is None or emit_completion_message is None:
        default_runner, default_emit = resolve_activity_pipeline_hooks()
        runner = runner or default_runner
        emit_completion_message = emit_completion_message or default_emit

    _sync_pipeline_logger(logger)

    return _run_activity_pipeline(
        cfg,
        options,
        runner=runner,
        emit_completion_message=emit_completion_message,
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run scripts.get_activity_data.

    Parameters
    ----------
    argv:
        Optional sequence of command-line arguments.

    Returns
    -------
    int
        Exit code returned by the script.
    """

    return _run("get_activity_data", argv)


__all__ = [
    "ActivityCommandOptions",
    "MIN_ACTIVITY_TIMEOUT",
    "run_activity_pipeline",
    "main",
]
