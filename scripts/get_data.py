"""Compatibility wrapper exposing :mod:`library.cli.commands.get_data`."""

from __future__ import annotations

import sys
from importlib import import_module

_module = import_module("library.cli.commands.get_data")
sys.modules[__name__] = _module

if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_module.main())
