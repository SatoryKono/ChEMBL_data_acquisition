"""Backwards compatible exports for the activity CLI command helpers."""

from __future__ import annotations

import argparse
import sys
from importlib import import_module
from types import ModuleType
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, Callable, Sequence

from library.cli.activity_api import (
    ActivityCommandOptions,
    MAX_ACTIVITY_CHUNK_SIZE,
    MIN_ACTIVITY_TIMEOUT,
    ensure_entrypoint_exports,
    main,
    register_activity_pipeline_hooks,
    resolve_activity_pipeline_hooks,
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

    try:
        entrypoint = _load_activity_entrypoint()
    except Exception:  # pragma: no cover - defensive fallback
        return _run_pipeline_fallback(*args, **kwargs)

    entrypoint_run = getattr(entrypoint, "run_chembl", None)
    if callable(entrypoint_run):
        return entrypoint_run(*args, **kwargs)

    return run_activity_pipeline(*args, **kwargs)


def _run_pipeline_fallback(*args: Any, **kwargs: Any) -> int:
    """Execute the activity pipeline without relying on the entrypoint module."""

    if not args:
        msg = "run_chembl requires a configuration object"
        raise TypeError(msg)

    cfg, *remaining = args
    if not remaining and "options" in kwargs:
        remaining = [kwargs.pop("options")]

    if not remaining:
        msg = "run_chembl requires activity options when falling back"
        raise TypeError(msg)

    options_source = remaining[0]
    options = _coerce_activity_options(options_source)

    runner, emit_completion = _resolve_pipeline_hooks_for_fallback()
    return run_activity_pipeline(
        cfg,
        options,
        runner=runner,
        emit_completion_message=emit_completion,
    )


def _coerce_activity_options(value: Any) -> ActivityCommandOptions:
    """Convert legacy CLI arguments into :class:`ActivityCommandOptions`."""

    if isinstance(value, ActivityCommandOptions):
        return value
    if isinstance(value, argparse.Namespace):
        return ActivityCommandOptions(
            input_csv=value.input_csv,
            output_csv=_namespace_value(value, "output_csv"),
            final_output=_coerce_final_output(value),
            limit=_namespace_value(value, "limit"),
            offset=int(_namespace_value(value, "offset", 0) or 0),
            timeout=_namespace_value(value, "timeout"),
            batch_size=_namespace_value(value, "batch_size"),
            workers=_namespace_value(value, "workers"),
            dry_run=bool(_namespace_value(value, "dry_run", False)),
            skip_existing=bool(_namespace_value(value, "skip_existing", False)),
            force=bool(_namespace_value(value, "force", False)),
            invocation=_coerce_invocation(value),
        )
    if isinstance(value, Mapping):
        return ActivityCommandOptions(**value)

    msg = "Unsupported activity options type for run_chembl fallback"
    raise TypeError(msg)


def _namespace_value(namespace: argparse.Namespace, name: str, default: Any | None = None) -> Any | None:
    """Return ``name`` from ``namespace`` ignoring ``argparse.SUPPRESS`` sentinels."""

    candidate = getattr(namespace, name, default)
    if candidate is argparse.SUPPRESS:
        return default
    return candidate


def _coerce_final_output(namespace: argparse.Namespace) -> Any | None:
    """Derive the final output path respecting legacy namespace attributes."""

    direct = _namespace_value(namespace, "final_output")
    if direct is not None:
        return direct
    legacy = _namespace_value(namespace, "final_out")
    if legacy is not None:
        return legacy
    return _namespace_value(namespace, "output_csv")


def _coerce_invocation(namespace: argparse.Namespace) -> Sequence[str] | None:
    """Normalise the ``invocation`` attribute from CLI namespaces."""

    invocation = _namespace_value(namespace, "invocation")
    if invocation is None:
        return None
    if isinstance(invocation, Sequence):
        return tuple(str(part) for part in invocation)
    return (str(invocation),)


def _resolve_pipeline_hooks_for_fallback() -> tuple[
    Callable[["Config", argparse.Namespace], int],
    Callable[..., None],
]:
    """Obtain pipeline hooks ensuring the fallback does not recurse."""

    runner, emit_completion = resolve_activity_pipeline_hooks()
    if runner is not run_chembl:
        return runner, emit_completion

    try:  # pragma: no cover - defensive guard for alternative packaging layouts
        from scripts import get_activity_data as legacy_cli  # type: ignore
    except Exception:  # pragma: no cover - maintain best-effort compatibility
        msg = "run_chembl fallback could not locate a pipeline runner"
        raise RuntimeError(msg)

    legacy_runner = getattr(legacy_cli, "run_chembl", None)
    legacy_emit = getattr(legacy_cli, "_emit_completion_message", None)
    if not callable(legacy_runner) or not callable(legacy_emit) or legacy_runner is run_chembl:
        msg = "run_chembl fallback could not locate a valid pipeline runner"
        raise RuntimeError(msg)

    register_activity_pipeline_hooks(
        runner=legacy_runner,
        emit_completion_message=legacy_emit,
    )
    return legacy_runner, legacy_emit


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
