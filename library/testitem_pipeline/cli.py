"""Legacy compatibility wrapper for test item pipeline CLI helpers.

Changelog
========
* Re-export implementation from :mod:`library.pipelines.testitem.cli` and warn
  about the deprecated import path.
"""

from __future__ import annotations

from warnings import warn

from library.pipelines.testitem.cli import *  # noqa: F401,F403

warn(
    "library.testitem_pipeline.cli is deprecated; import from "
    "library.pipelines.testitem.cli instead.",
    DeprecationWarning,
    stacklevel=2,
)

try:
    from library.pipelines.testitem.cli import __all__ as _CLI_ALL
except ImportError:  # pragma: no cover - defensive fallback
    __all__ = [name for name in globals() if not name.startswith("_")]
else:
    __all__ = list(_CLI_ALL)
