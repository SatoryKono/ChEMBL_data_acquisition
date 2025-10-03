"""Legacy compatibility wrapper for test item catalog helpers.

Changelog
========
* Re-export implementation from :mod:`library.pipelines.testitem.catalog` and
  warn about the deprecated import path.
"""

from __future__ import annotations

from warnings import warn

from library.pipelines.testitem.catalog import *  # noqa: F401,F403
from library.pipelines.testitem.catalog import __all__ as _CATALOG_ALL

warn(
    "library.testitem_pipeline.catalog is deprecated; use "
    "library.pipelines.testitem.catalog instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = list(_CATALOG_ALL)
