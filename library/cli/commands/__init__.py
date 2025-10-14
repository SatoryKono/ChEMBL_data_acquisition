"""Shared helpers for CLI command wrappers."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from importlib import import_module
from inspect import signature
from types import ModuleType
from typing import Final

CommandName = str

_COMMAND_MODULE_PATHS: Final[dict[CommandName, str]] = {
    "check_determinism": "library.utils.cli_tools.check_determinism",
    "chunk_io": "library.utils.cli_tools.chunk_io",
    "csv_utils": "library.utils.cli_tools.csv_utils",
    "get_activities": "library.utils.cli_tools.get_activities",
    "get_activity_data": "library.cli.commands.get_activity_data",
    "get_assay_data": "library.cli.commands.get_assay_data",
    "get_data": "library.cli.commands.get_data",
    "get_document_data": "library.cli.commands.get_document_data",
    "get_document_type": "library.utils.cli_tools.get_document_type",
    "get_input_initialisation": "library.utils.cli_tools.get_input_initialisation",
    "get_target_data": "library.cli.commands.get_target_data",
    "get_testitem_data": "library.cli.commands.get_testitem_data",
    "mapper": "library.utils.cli_tools.mapper",
    "table_quality": "library.utils.cli_tools.table_quality",
}


def resolve_command_path(module: CommandName) -> str:
    """Return the fully-qualified module path for ``module``."""

    if "." in module:
        return module

    try:
        return _COMMAND_MODULE_PATHS[module]
    except KeyError as exc:  # pragma: no cover - defensive guard
        msg = f"Could not resolve CLI module '{module}'"
        raise ModuleNotFoundError(msg) from exc


def _resolve_module(module: str) -> ModuleType:
    """Return the module object for the requested CLI tool."""

    module_path = resolve_command_path(module)
    try:
        return import_module(module_path)
    except ModuleNotFoundError as exc:  # pragma: no cover - defensive branch
        if module_path == module:
            raise
        msg = f"Could not import '{module_path}' for CLI module '{module}'"
        raise ModuleNotFoundError(msg) from exc


def _run(module: str, argv: Sequence[str] | None = None) -> int:
    """Execute ``module.main`` from the supported CLI namespaces."""

    main_func: Callable[..., int] = _resolve_module(module).main
    params = signature(main_func).parameters
    if not params:
        return main_func()
    if argv is None:
        argv = []
    return main_func(argv)


__all__: list[str] = ["_run", "resolve_command_path"]
