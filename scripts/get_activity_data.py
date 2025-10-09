"""Thin wrapper exposing the activity pipeline CLI entry point."""

from __future__ import annotations

if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

from importlib import import_module
import sys

from library.cli.entrypoints.activity import (
    ActivityPipelineCLI,
    DEFAULT_INPUT_NAME,
    DEFAULT_OUTPUT_STEM,
    PROGRAM_NAME,
    _emit_completion_message,
    build_parser,
    main,
    run,
    run_chembl,
)
from library.pipelines.activity import runner as _activity_runner
from library.pipelines.activity.runner import (
    MAX_ACTIVITY_CHUNK_SIZE,
    register_activity_pipeline_hooks,
)

_activity = import_module("library.cli.entrypoints.activity")

register_activity_pipeline_hooks(
    runner=run_chembl,
    emit_completion_message=_emit_completion_message,
)

setattr(_activity, "MAX_ACTIVITY_CHUNK_SIZE", MAX_ACTIVITY_CHUNK_SIZE)


class _LoggerProxy:
    """Proxy ``library.pipelines.activity.runner`` logger to the CLI module."""

    def __getattr__(self, name: str):  # type: ignore[override]
        return getattr(_activity, "logger").__getattribute__(name)


_activity_runner.logger = _LoggerProxy()

_extended_exports = (
    "ActivityPipelineCLI",
    "DEFAULT_INPUT_NAME",
    "DEFAULT_OUTPUT_STEM",
    "PROGRAM_NAME",
    "run",
    "run_chembl",
    "_emit_completion_message",
    "main",
    "build_parser",
    "MAX_ACTIVITY_CHUNK_SIZE",
)

_activity.__all__ = tuple(
    dict.fromkeys(getattr(_activity, "__all__", tuple()) + _extended_exports)
)

sys.modules[__name__] = _activity
sys.modules.setdefault("scripts.get_activity_data", _activity)

if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_activity.main())
