"""Compatibility wrapper for :mod:`library.testitem_pipeline.cli`."""

from __future__ import annotations

import warnings
from importlib import import_module

_modern_cli = import_module("library.pipelines.testitem.cli")

warnings.warn(
    "library.testitem_pipeline.cli is deprecated; use library.pipelines.testitem.cli instead.",
    DeprecationWarning,
    stacklevel=2,
)

_SKIP_COPY = {"__doc__", "__loader__", "__name__", "__package__", "__spec__", "__file__"}

for _name in dir(_modern_cli):
    if _name in _SKIP_COPY:
        continue
    globals()[_name] = getattr(_modern_cli, _name)

if "augment_pubchem" not in globals():
    from library.pipelines.testitem.pubchem import augment_pubchem as _augment

    globals()["augment_pubchem"] = _augment

__all__ = [name for name in globals() if not name.startswith("__")]
