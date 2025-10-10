"""Domain-specific data pipelines grouped by entity type.

Changelog
========
* Replace eager imports with lazy module loaders to avoid circular
  dependencies during initialisation.
"""

from __future__ import annotations

# ===== Modules =====
from importlib import import_module
from typing import Any

# ===== Helpers =====
_LAZY_SUBMODULES = {
    "activity",
    "assay",
    "cellline",
    "common",
    "document",
    "target",
    "testitem",
    "tissue",
}


def _load_submodule(name: str) -> Any:
    """Import ``library.pipelines.<name>`` and cache the module in globals."""

    module = import_module(f"{__name__}.{name}")
    globals()[name] = module
    return module


# ===== Exports =====
__all__ = [
    "activity",
    "assay",
    "cellline",
    "common",
    "document",
    "target",
    "testitem",
    "tissue",
]


def __getattr__(name: str) -> Any:
    """Expose pipeline submodules on demand to break import cycles."""

    if name in _LAZY_SUBMODULES:
        return _load_submodule(name)
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def __dir__() -> list[str]:
    """Provide completion support for the exported submodules."""

    return sorted({*globals(), *__all__, *_LAZY_SUBMODULES})
