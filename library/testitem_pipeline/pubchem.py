"""Compatibility wrapper for :mod:`library.testitem_pipeline.pubchem`."""

from __future__ import annotations

import warnings
from importlib import import_module

_modern_pubchem = import_module("library.pipelines.testitem.pubchem")

warnings.warn(
    "library.testitem_pipeline.pubchem is deprecated; use library.pipelines.testitem.pubchem instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = list(getattr(_modern_pubchem, "__all__", ()))

globals().update({name: getattr(_modern_pubchem, name) for name in __all__})
