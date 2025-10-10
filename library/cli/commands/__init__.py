"""Shared helpers for CLI command wrappers."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from importlib import import_module
from inspect import signature
from types import ModuleType


def _resolve_module(module: str) -> ModuleType:
    """Return the module object for the requested CLI tool."""

    if "." in module:
        return import_module(module)

    last_error: ModuleNotFoundError | None = None
    for candidate in (
        "library.cli.commands." + module,
        "library.utils.cli_tools." + module,
        "scripts." + module,
    ):
        try:
            return import_module(candidate)
        except ModuleNotFoundError as exc:  # pragma: no cover - defensive branch
            last_error = exc

    if last_error is None:  # pragma: no cover - should never happen
        msg = f"Could not resolve CLI module '{module}'"
        raise ModuleNotFoundError(msg)
    raise last_error


def _run(module: str, argv: Sequence[str] | None = None) -> int:
    """Execute ``module.main`` from the supported CLI namespaces."""

    main_func: Callable[..., int] = _resolve_module(module).main
    params = signature(main_func).parameters
    if not params:
        return main_func()
    if argv is None:
        argv = []
    return main_func(argv)


__all__: list[str] = ["_run"]
