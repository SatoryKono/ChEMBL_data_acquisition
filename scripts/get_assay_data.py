"""Compatibility wrapper exposing the assay pipeline CLI entry point."""

from __future__ import annotations

import sys
from importlib import import_module
from types import ModuleType

# ruff: noqa: E402 - bootstrap adjusts import order for script execution
if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a package module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli


def _load_assay_module() -> ModuleType:
    """Return the authoritative assay CLI module."""

    return import_module("library.cli.entrypoints.assay")


_module = _load_assay_module()
__all__ = tuple(getattr(_module, "__all__", ()))


def __getattr__(name: str) -> object:
    """Delegate attribute access to :mod:`library.cli.entrypoints.assay`."""

    try:
        return getattr(_module, name)
    except AttributeError as exc:  # pragma: no cover - propagate missing attrs
        raise AttributeError(name) from exc


sys.modules[__name__] = _module
sys.modules.setdefault("scripts.get_assay_data", _module)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_module.main())
