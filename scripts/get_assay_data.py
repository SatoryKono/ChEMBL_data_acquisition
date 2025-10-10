"""Compatibility wrapper exposing :mod:`library.cli.commands.get_assay_data`."""

from __future__ import annotations

import sys
from importlib import import_module
from typing import TYPE_CHECKING

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if TYPE_CHECKING:
    from . import _bootstrap as _bootstrap_module
elif __package__ in {None, ""}:
    import _bootstrap as _bootstrap_module  # pragma: no cover - CLI fallback
else:  # pragma: no cover - executed when imported as a module
    from . import _bootstrap as _bootstrap_module

bootstrap_cli = _bootstrap_module.bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli
del _bootstrap_module

_MODULE = import_module("library.cli.commands.get_assay_data")
sys.modules[__name__] = _MODULE

if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_MODULE.main())
