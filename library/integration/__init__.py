"""Adapters and façade modules wrapping third-party integrations.

Changelog
========
* Replace eager submodule imports with lazy loading to avoid circular
  dependencies when accessed during the test item pipeline bootstrap.
"""

from __future__ import annotations

# ===== Modules =====
from importlib import import_module
from typing import Any

# ===== Helpers =====
_LAZY_SUBMODULES = {
    "input_initialisation_library",
    "mapper_batch_library",
    "mapper_library",
    "uniprot_library",
}


def _load_submodule(name: str) -> Any:
    """Import ``library.integration.<name>`` and memoise the module object."""

    module = import_module(f"{__name__}.{name}")
    globals()[name] = module
    return module


# ===== Exports =====
__all__ = [
    "input_initialisation_library",
    "mapper_batch_library",
    "mapper_library",
    "uniprot_library",
]


def __getattr__(name: str) -> Any:
    """Resolve lazy integration submodules on first attribute access."""

    if name in _LAZY_SUBMODULES:
        return _load_submodule(name)
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def __dir__() -> list[str]:
    """Expose eager and lazy exports for introspection tools."""

    return sorted({*globals(), *__all__, *_LAZY_SUBMODULES})
