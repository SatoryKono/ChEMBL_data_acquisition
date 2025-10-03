"""Compatibility wrapper for :mod:`library.testitem_pipeline.catalog`."""

from __future__ import annotations

import warnings
from importlib import import_module

_modern_catalog = import_module("library.pipelines.testitem.catalog")

warnings.warn(
    "library.testitem_pipeline.catalog is deprecated; use library.pipelines.testitem.catalog instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = list(getattr(_modern_catalog, "__all__", ()))

globals().update({name: getattr(_modern_catalog, name) for name in __all__})
