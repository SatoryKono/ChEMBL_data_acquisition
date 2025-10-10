"""Backwards compatible exports for the activity CLI command helpers."""

from __future__ import annotations

import argparse
import sys
from importlib import import_module
from types import ModuleType
from typing import TYPE_CHECKING, Any

from library.cli.activity_api import (
    ActivityCommandOptions,
    MAX_ACTIVITY_CHUNK_SIZE,
    MIN_ACTIVITY_TIMEOUT,
    ensure_entrypoint_exports,
    main,
    register_activity_pipeline_hooks,
    run_activity_pipeline,
)


_BASE_EXPORTS = (
    "ActivityCommandOptions",
    "MAX_ACTIVITY_CHUNK_SIZE",
    "MIN_ACTIVITY_TIMEOUT",
    "register_activity_pipeline_hooks",
    "run_activity_pipeline",
    "main",
    "run_chembl",
    "_emit_completion_message",
)
__all__ = _BASE_EXPORTS


def _bootstrap_entrypoint_exports() -> None:
    """Ensure the entrypoint module exposes the shared command helpers."""

    try:
        ensure_entrypoint_exports()
    except Exception:  # pragma: no cover - defensive guard
        # Import side effects should never break command module import.
        pass


_bootstrap_entrypoint_exports()


def _synchronise_wrapper_module() -> None:
    """Ensure :mod:`scripts.get_activity_data` aliases this module when loaded."""

    wrapper_name = "scripts.get_activity_data"
    self_module = sys.modules[__name__]
    wrapper_module = sys.modules.get(wrapper_name)
    if wrapper_module is None or wrapper_module is self_module:
        return
    if getattr(wrapper_module, "_MODULE", None) is self_module:
        sys.modules[wrapper_name] = self_module


_synchronise_wrapper_module()


def _load_activity_entrypoint() -> ModuleType:
    """Return the activity entrypoint module ensuring compatibility exports."""

    module = import_module("library.cli.entrypoints.activity")
    ensure_entrypoint_exports(module)
    return module


if TYPE_CHECKING:  # pragma: no cover - imported for type checking only
    from library.config import Config


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the legacy :mod:`scripts.get_activity_data` runner."""

    entrypoint = _load_activity_entrypoint()
    return entrypoint.run_chembl(cfg, args)


def _emit_completion_message(*args: Any, **kwargs: Any) -> None:
    """Forward completion events to the entrypoint implementation."""

    entrypoint = _load_activity_entrypoint()
    entrypoint._emit_completion_message(*args, **kwargs)


def __getattr__(name: str) -> Any:
    """Defer missing attributes to :mod:`library.cli.entrypoints.activity`."""

    entrypoint = _load_activity_entrypoint()
    try:
        return getattr(entrypoint, name)
    except AttributeError as exc:  # pragma: no cover - passthrough helper
        raise AttributeError(name) from exc


def __dir__() -> list[str]:  # pragma: no cover - helper for introspection
    entrypoint = _load_activity_entrypoint()
    return sorted({*globals().keys(), *dir(entrypoint)})


def _entrypoint_public_names() -> tuple[str, ...]:
    """Return the exported public names from the entrypoint module."""

    try:
        entrypoint = _load_activity_entrypoint()
    except Exception:  # pragma: no cover - defensive fallback
        return ()
    exported = getattr(entrypoint, "__all__", ())
    return tuple(str(name) for name in exported)


__all__ = tuple(dict.fromkeys((*_BASE_EXPORTS, *_entrypoint_public_names())))
