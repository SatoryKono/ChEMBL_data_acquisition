"""Deprecated catalog helpers for :mod:`library.pipelines.testitem`."""

from __future__ import annotations

import warnings

from library.pipelines.testitem import catalog as _catalog
from library.pipelines.testitem.core import _CATALOG_EXPORTS

warnings.warn(
    "library.testitem_pipeline.catalog is deprecated; import from "
    "library.pipelines.testitem.catalog instead.",
    DeprecationWarning,
    stacklevel=2,
)

for _name in _CATALOG_EXPORTS:
    globals()[_name] = getattr(_catalog, _name)
del _name

__all__ = list(_CATALOG_EXPORTS)
