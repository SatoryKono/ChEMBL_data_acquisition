"""CLI entry point implementations for pipeline scripts."""

from __future__ import annotations

from importlib import import_module
from typing import Any, TYPE_CHECKING

__all__ = ["ActivityPipelineCLI"]


if TYPE_CHECKING:  # pragma: no cover - imported only for type checking
    from .activity import ActivityPipelineCLI as ActivityPipelineCLITyped


def __getattr__(name: str) -> Any:
    """Lazily import heavy entry points on demand."""

    if name == "ActivityPipelineCLI":
        module = import_module("library.cli.entrypoints.activity")
        return getattr(module, name)
    msg = f"module 'library.cli.entrypoints' has no attribute {name!r}"
    raise AttributeError(msg)
