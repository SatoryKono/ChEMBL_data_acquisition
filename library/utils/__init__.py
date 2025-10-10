"""Utility helpers shared across the ChEMBL acquisition tooling."""

from __future__ import annotations

from importlib import import_module
from typing import Any

from . import atomic, bootstrap
from .data_correlation import build_correlation_matrix
from .qc_report import build_qc_summary

__all__ = [
    "atomic",
    "bootstrap",
    "config",
    "build_qc_summary",
    "build_correlation_matrix",
]


def __getattr__(name: str) -> Any:
    """Lazily import heavy utility submodules."""

    if name == "config":
        module = import_module(f"{__name__}.config")
        globals()[name] = module
        return module
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def __dir__() -> list[str]:
    return sorted({*globals(), *__all__})
