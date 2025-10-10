from __future__ import annotations

import argparse
import importlib
import sys
from functools import wraps
from collections.abc import Callable, Sequence
from types import ModuleType
from typing import cast

from library.common.log import logger
from library.common.logging_setup import Logger
from library.config import Config
from library.pipelines.activity.runner import (
    MAX_ACTIVITY_CHUNK_SIZE,
    MIN_ACTIVITY_TIMEOUT,
    ActivityCommandOptions,
    register_activity_pipeline_hooks,
    resolve_activity_pipeline_hooks,
    run_activity_pipeline as _run_activity_pipeline,
)

def _sync_pipeline_logger(current_logger: Logger) -> None:
    """Align the shared pipeline logger with ``current_logger`` when possible."""

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
    """Delegate to :mod:`library.cli.entrypoints.activity` without circular imports."""

    activity_entrypoint = importlib.import_module(
        "library.cli.entrypoints.activity"
    )
    entry_main = cast(Callable[[Sequence[str] | None], int], getattr(activity_entrypoint, "main"))
    return entry_main(argv)


_ENTRYPOINT_MODULE_NAME = "library.cli.entrypoints.activity"
_ENTRYPOINT_MODULE_CACHE: ModuleType | None = None


def _load_entrypoint_module() -> ModuleType:
    global _ENTRYPOINT_MODULE_CACHE

    if _ENTRYPOINT_MODULE_CACHE is None:
        _ENTRYPOINT_MODULE_CACHE = importlib.import_module(_ENTRYPOINT_MODULE_NAME)
    return _ENTRYPOINT_MODULE_CACHE


_COMPAT_EXPORTS: tuple[str, ...] = (
    "ActivityCommandOptions",
    "MAX_ACTIVITY_CHUNK_SIZE",
    "MIN_ACTIVITY_TIMEOUT",
    "register_activity_pipeline_hooks",
    "run_activity_pipeline",
)

__all__ = [*_COMPAT_EXPORTS, "main"]


def _publish_compat_exports() -> ModuleType:
    entrypoint = _load_entrypoint_module()

    for name in _COMPAT_EXPORTS:
        setattr(entrypoint, name, globals()[name])

    entrypoint_all = getattr(entrypoint, "__all__", ())
    if entrypoint_all:
        combined = tuple(dict.fromkeys((*entrypoint_all, *_COMPAT_EXPORTS)))
        entrypoint.__all__ = combined  # type: ignore[attr-defined]

    original_emit = getattr(entrypoint, "_emit_completion_message", None)
    if original_emit is not None and not getattr(original_emit, "_commands_skip_wrapper", False):

        @wraps(original_emit)
        def _emit_completion_message_with_skip(*args, **kwargs):
            result = original_emit(*args, **kwargs)

            mode = kwargs.get("mode")
            output_path = kwargs.get("output_path")
            if mode == "skip_existing":
                payload = {
                    "output": str(output_path) if output_path is not None else None,
                    "rows": int(kwargs.get("processed_rows") or 0),
                    "duration_s": float(kwargs.get("duration_s", 0.0)),
                    "mode": mode,
                }
                logger_obj = getattr(entrypoint, "logger", None)
                if logger_obj is not None:
                    try:
                        logger_obj.info("activity_pipeline_completion", **payload)
                    except Exception:  # pragma: no cover - defensive guard for non-structlog loggers
                        pass

            return result

        setattr(_emit_completion_message_with_skip, "_commands_skip_wrapper", True)
        setattr(entrypoint, "_emit_completion_message", _emit_completion_message_with_skip)

    setattr(entrypoint, "_activity_cli_commands", entrypoint)
    return entrypoint


_ENTRYPOINT_PROXY = _publish_compat_exports()
sys.modules[__name__] = _ENTRYPOINT_PROXY
