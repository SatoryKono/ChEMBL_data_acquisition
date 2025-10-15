"""Compatibility wrappers for duplicate analysis tooling.

The quality assurance utilities moved to :mod:`library.qa` when the project
was reorganised.  This module re-exports the public entry points so that
imports such as ``import project.duplicate_analysis.table_quality`` keep
working until dependent projects migrate.
"""

from __future__ import annotations

from importlib import import_module
from types import ModuleType
from typing import Final

_REDIRECTS: Final[dict[str, str]] = {
    "table_quality": "library.qa.table_quality",
    "reporting": "library.qa.reporting",
    "validation": "library.qa.validation",
}

_cache: dict[str, ModuleType] = {}


def _load(name: str) -> ModuleType:
    try:
        target = _REDIRECTS[name]
    except KeyError as exc:  # pragma: no cover - defensive programming
        raise AttributeError(
            f"module 'project.duplicate_analysis' has no attribute '{name}'"
        ) from exc
    module = import_module(target)
    _cache[name] = module
    return module


def __getattr__(name: str) -> ModuleType:
    return _cache.get(name) or _load(name)


def __dir__() -> list[str]:  # pragma: no cover - reflective helper
    return sorted(set(globals()) | set(_REDIRECTS))


__all__ = ["_load"]
