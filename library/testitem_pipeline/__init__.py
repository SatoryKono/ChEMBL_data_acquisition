"""Compatibility layer for the historic :mod:`library.testitem_pipeline` namespace."""

from __future__ import annotations

import warnings
from importlib import import_module

_modern = import_module("library.pipelines.testitem")

warnings.warn(
    "library.testitem_pipeline is deprecated; use library.pipelines.testitem instead.",
    DeprecationWarning,
    stacklevel=2,
)

_SKIP_COPY = {"__doc__", "__loader__", "__name__", "__package__", "__spec__", "__file__", "__path__"}

for _name in dir(_modern):
    if _name in _SKIP_COPY:
        continue
    globals()[_name] = getattr(_modern, _name)

__all__ = list(getattr(_modern, "__all__", ()))
