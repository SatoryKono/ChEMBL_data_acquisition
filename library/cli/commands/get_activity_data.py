"""Backwards compatible exports for the activity CLI command helpers."""

from __future__ import annotations

from library.cli.activity_api import (
    ActivityCommandOptions,
    MAX_ACTIVITY_CHUNK_SIZE,
    MIN_ACTIVITY_TIMEOUT,
    ensure_entrypoint_exports,
    main,
    register_activity_pipeline_hooks,
    run_activity_pipeline,
)


__all__ = (
    "ActivityCommandOptions",
    "MAX_ACTIVITY_CHUNK_SIZE",
    "MIN_ACTIVITY_TIMEOUT",
    "register_activity_pipeline_hooks",
    "run_chembl",
    "run_activity_pipeline",
    "main",
)


def run_chembl(*args, **kwargs):
    """Backward compatible alias delegating to :func:`run_activity_pipeline`."""

    return run_activity_pipeline(*args, **kwargs)


def _bootstrap_entrypoint_exports() -> None:
    """Ensure the entrypoint module exposes the shared command helpers."""

    try:
        ensure_entrypoint_exports()
    except Exception:  # pragma: no cover - defensive guard
        # Import side effects should never break command module import.
        pass


_bootstrap_entrypoint_exports()
