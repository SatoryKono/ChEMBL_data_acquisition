"""Helpers for resolving dotted import paths for pipeline configuration."""
from __future__ import annotations

import importlib
from typing import Any

from .types import ImportResolutionError


def resolve_dotted_path(path: str) -> Any:
    """Return the object referenced by ``path``.

    Parameters
    ----------
    path:
        Dotted path in the form ``"package.module:attribute"`` or
        ``"package.module.attribute"``.
    """

    module_path: str
    attribute: str | None

    if ":" in path:
        module_path, attribute = path.split(":", 1)
    else:
        parts = path.split(".")
        if len(parts) < 2:
            raise ImportResolutionError(f"Invalid dotted path '{path}'")
        module_path = ".".join(parts[:-1])
        attribute = parts[-1]

    try:
        module = importlib.import_module(module_path)
    except ModuleNotFoundError as exc:  # pragma: no cover - defensive
        raise ImportResolutionError(f"Cannot import module '{module_path}'") from exc

    if attribute is None:
        return module

    try:
        return getattr(module, attribute)
    except AttributeError as exc:  # pragma: no cover - defensive
        raise ImportResolutionError(
            f"Module '{module_path}' has no attribute '{attribute}'"
        ) from exc


__all__ = ["resolve_dotted_path"]
