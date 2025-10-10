"""Compatibility wrapper exposing :mod:`library.cli.commands.get_data`."""

from __future__ import annotations

import sys
from importlib import import_module

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

_module = import_module("library.cli.commands.get_data")
sys.modules[__name__] = _module
sys.modules["scripts.get_data"] = _module

if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_module.main())
