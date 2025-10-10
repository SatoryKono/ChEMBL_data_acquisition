"""Compatibility wrapper exposing :mod:`library.cli.commands.get_assay_data`."""

from __future__ import annotations

import argparse
import sys
from importlib import import_module

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

_module = import_module("library.cli.commands.get_assay_data")

__all__ = getattr(_module, "__all__", [])

_COMPAT_EXPORTS = {
    "ASSAY_COLUMNS",
    "ASSAY_OUTPUT_DROP_COLUMNS",
    "_ASSAY_OUTPUT_DROP_COLUMNS",
    "ASSAY_MAX_IDS_PER_REQUEST",
    "_ASSAY_MAX_IDS_PER_REQUEST",
    "DEFAULT_INPUT_NAME",
    "DEFAULT_OUTPUT_STEM",
    "MAX_ASSAY_CHUNK_SIZE",
}

for _name in set(__all__) | _COMPAT_EXPORTS:
    if hasattr(_module, _name):
        globals()[_name] = getattr(_module, _name)


def __getattr__(name: str):  # pragma: no cover - passthrough helper
    return getattr(_module, name)


def __dir__():  # pragma: no cover - passthrough helper
    return sorted(set(globals()) | set(dir(_module)))


sys.modules[__name__] = _module
sys.modules.setdefault("scripts.get_assay_data", _module)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_module.main())
