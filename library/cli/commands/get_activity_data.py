from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence

from library.config import Config
from library.common.log import logger
from library.pipelines.activity.runner import (
    ActivityCommandOptions,
    MIN_ACTIVITY_TIMEOUT,
    resolve_activity_pipeline_hooks,
    run_activity_pipeline as _run_activity_pipeline,
)

from . import _run

def _sync_pipeline_logger(logger: object) -> None:
    """Align the shared pipeline logger with ``logger`` when possible."""

    try:
        from library.common import log as _common_log

        _common_log.logger = logger  # type: ignore[assignment]
    except Exception:  # pragma: no cover - defensive
        pass

    try:
        from library.pipelines.activity import runner as activity_runner

        activity_runner.logger = logger  # type: ignore[assignment]
    except Exception:  # pragma: no cover - defensive
        pass


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
