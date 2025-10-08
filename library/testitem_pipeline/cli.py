"""Deprecated CLI helpers for :mod:`library.pipelines.testitem`."""

from __future__ import annotations

import warnings

from library.pipelines.testitem import cli as _cli
from library.pipelines.testitem.core import _CLI_EXPORTS

warnings.warn(
    "library.testitem_pipeline.cli is deprecated; import from "
    "library.pipelines.testitem.cli instead.",
    DeprecationWarning,
    stacklevel=2,
)

for _name in _CLI_EXPORTS:
    globals()[_name] = getattr(_cli, _name)
del _name

__all__ = list(_CLI_EXPORTS)
