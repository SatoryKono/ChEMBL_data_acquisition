"""Backwards compatible exports for the activity CLI command helpers."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable, Mapping
from importlib import import_module
from types import ModuleType
from typing import TYPE_CHECKING, Any

import library.cli.logging as cli_logging
from library.cli.activity_api import (
    ensure_entrypoint_exports,
    run_activity_pipeline,
)
from library.cli_utils import run_cli_command as _run_cli_command

_BASE_EXPORTS = (
    "ActivityCommandOptions",
    "MAX_ACTIVITY_CHUNK_SIZE",
    "MIN_ACTIVITY_TIMEOUT",
    "register_activity_pipeline_hooks",
    "run_activity_pipeline",
    "main",
    "run_cli_command",
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
    from library.cli import Logger, LoggerConfig
    from library.config import Config


def _normalise_date_token(value: object) -> str | None:
    """Return a stripped date token or ``None`` when not provided."""

    if not isinstance(value, str):
        return None
    token = value.strip()
    return token or None


def _resolve_configured_prefix(cfg: Config) -> str | None:
    """Return the configured default date prefix when available."""

    io_cfg = getattr(cfg, "io", None)
    if io_cfg is None:
        local_cfg = getattr(cfg, "local", None)
        io_cfg = getattr(local_cfg, "io", None)

    candidate = getattr(io_cfg, "default_date_prefix", None) if io_cfg else None
    if isinstance(candidate, str):
        candidate = candidate.strip() or None
    else:
        candidate = None
    return candidate


def _resolve_stamp_mode(cfg: Config, args: argparse.Namespace) -> str | None:
    """Return the effective output stamp mode for CLI invocation."""

    explicit = getattr(args, "output_stamp_mode", None)
    if isinstance(explicit, str) and explicit.strip():
        return explicit.strip()

    io_cfg = getattr(cfg, "io", None)
    if io_cfg is None:
        local_cfg = getattr(cfg, "local", None)
        io_cfg = getattr(local_cfg, "io", None)

    stamp_mode = getattr(io_cfg, "output_stamp_mode", None) if io_cfg else None
    if isinstance(stamp_mode, str):
        return stamp_mode.strip() or None
    return None


def _ensure_default_date(cfg: Config, args: argparse.Namespace) -> str | None:
    """Populate ``args.date`` honouring configuration defaults."""

    explicit = _normalise_date_token(getattr(args, "date", None))
    if explicit is not None:
        args.date = explicit
        return explicit

    configured = _resolve_configured_prefix(cfg)
    if configured is not None:
        args.date = configured
        return configured

    stamp_mode = _resolve_stamp_mode(cfg, args)
    if stamp_mode == "require":
        fallback = cli_logging._current_date_str()
        args.date = fallback
        return fallback

    if hasattr(args, "date"):
        args.date = None
    return None


def run_chembl(*args: Any, **kwargs: Any) -> int:
    """Execute the activity pipeline using the entrypoint when available.

    Falls back to :func:`run_activity_pipeline` when the entrypoint cannot be
    imported, preserving backwards compatibility with older deployments that
    relied on the direct pipeline invocation helper.
    """

    try:
        entrypoint = _load_activity_entrypoint()
    except Exception:  # pragma: no cover - defensive fallback
        return run_activity_pipeline(*args, **kwargs)

    entrypoint_run = getattr(entrypoint, "run_chembl", None)
    if callable(entrypoint_run):
        return entrypoint_run(*args, **kwargs)

    return run_activity_pipeline(*args, **kwargs)


def _emit_completion_message(*args: Any, **kwargs: Any) -> None:
    """Forward completion events to the entrypoint implementation."""

    entrypoint = _load_activity_entrypoint()
    entrypoint._emit_completion_message(*args, **kwargs)


def run_cli_command(
    *,
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
    base_parser: argparse.ArgumentParser | None = None,
    log_cfg: LoggerConfig,
    mapping: Mapping[str, str],
    run: Callable[[Config, argparse.Namespace], int],
    logger: Logger | None = None,
    postprocess_enabled: bool | None = None,
) -> int:
    """Wrap :func:`library.cli_utils.run_cli_command` adding date defaults."""

    def _run_with_date(cfg: Config, parsed_args: argparse.Namespace) -> int:
        _ensure_default_date(cfg, parsed_args)
        return run(cfg, parsed_args)

    return _run_cli_command(
        args=args,
        parser=parser,
        base_parser=base_parser,
        log_cfg=log_cfg,
        mapping=mapping,
        run=_run_with_date,
        logger=logger,
        postprocess_enabled=postprocess_enabled,
    )


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
