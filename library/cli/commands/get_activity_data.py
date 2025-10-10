"""Helpers exposing the activity pipeline CLI under :mod:`library.cli.commands`."""
from __future__ import annotations

import argparse
import importlib
from collections.abc import Callable, Sequence
from types import ModuleType
from typing import cast

from library.common.log import logger
from library.config import Config
from library.pipelines.activity.runner import (
    MIN_ACTIVITY_TIMEOUT,
    ActivityCommandOptions,
    resolve_activity_pipeline_hooks,
)
from library.pipelines.activity.runner import (
    run_activity_pipeline as _run_activity_pipeline,
)


def _sync_pipeline_logger(current_logger: object) -> None:
    """Align shared pipeline loggers with ``current_logger`` when possible."""

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
    """Invoke the primary activity CLI entry point."""

    from library.cli.entrypoints import activity as activity_cli

    return activity_cli.main(argv)


__all__ = [
    "ActivityCommandOptions",
    "MIN_ACTIVITY_TIMEOUT",
    "run_activity_pipeline",
    "main",
]
