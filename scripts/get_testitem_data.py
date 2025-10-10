"""Compatibility wrapper exposing :mod:`library.cli.commands.get_testitem_data`."""

from __future__ import annotations

import sys
from importlib import import_module

try:
    from _bootstrap import bootstrap_cli
except ImportError:  # pragma: no cover - namespace package fallback
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

_module = import_module("library.cli.commands.get_testitem_data")
sys.modules[__name__] = _module

if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_module.main())
