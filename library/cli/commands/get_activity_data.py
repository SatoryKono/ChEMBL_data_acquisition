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


def run_chembl(*args: Any, **kwargs: Any) -> int:
    """Execute the activity pipeline using the entrypoint when available.

    Falls back to :func:`run_activity_pipeline` when the entrypoint cannot be
    imported, preserving backwards compatibility with older deployments that
    relied on the direct pipeline invocation helper.
    """

    def _coerce_activity_options(value: Any) -> Any:
        if isinstance(value, ActivityCommandOptions):
            return value
        if isinstance(value, argparse.Namespace):
            output_csv = getattr(value, "output_csv", None)
            final_output = getattr(value, "final_output", None)
            if final_output is None:
                final_output = getattr(value, "final_out", None)
            return ActivityCommandOptions(
                input_csv=getattr(value, "input_csv"),
                output_csv=output_csv,
                final_output=final_output,
                limit=getattr(value, "limit", None),
                offset=getattr(value, "offset", 0),
                timeout=getattr(value, "timeout", None),
                batch_size=getattr(value, "batch_size", None),
                workers=getattr(value, "workers", None),
                dry_run=getattr(value, "dry_run", False),
                skip_existing=getattr(value, "skip_existing", False),
                force=getattr(value, "force", False),
                invocation=getattr(value, "invocation", None),
            )
        return value

    def _normalise_arguments(
        call_args: tuple[Any, ...], call_kwargs: dict[str, Any]
    ) -> tuple[tuple[Any, ...], dict[str, Any]]:
        if call_args:
            normalised_args = tuple(_coerce_activity_options(arg) for arg in call_args)
        else:
            normalised_args = call_args

        if call_kwargs:
            normalised_kwargs = {
                key: _coerce_activity_options(value)
                for key, value in call_kwargs.items()
            }
        else:
            normalised_kwargs = call_kwargs

        return normalised_args, normalised_kwargs

    try:
        entrypoint = _load_activity_entrypoint()
    except Exception:  # pragma: no cover - defensive fallback
        normalised_args, normalised_kwargs = _normalise_arguments(args, kwargs)
        return run_activity_pipeline(*normalised_args, **normalised_kwargs)

    entrypoint_run = getattr(entrypoint, "run_chembl", None)
    if callable(entrypoint_run):
        return entrypoint_run(*args, **kwargs)

    normalised_args, normalised_kwargs = _normalise_arguments(args, kwargs)
    return run_activity_pipeline(*normalised_args, **normalised_kwargs)


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
