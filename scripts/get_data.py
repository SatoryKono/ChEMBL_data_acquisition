"""Compatibility shim preserving the historical ``scripts.get_data`` module."""
from __future__ import annotations

from importlib import import_module
from types import ModuleType

import sys

# ruff: noqa: E402  # keep bootstrap import order for direct script execution
if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a package module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli


def _load_module() -> ModuleType:
    return import_module("library.cli.commands.get_data")


_module = _load_module()
__all__ = tuple(getattr(_module, "__all__", ()))

sys.modules[__name__] = _module
sys.modules.setdefault("scripts.get_data", _module)


def __getattr__(name: str) -> object:
    try:
        return getattr(_module, name)
    except AttributeError as exc:  # pragma: no cover - passthrough for missing attrs
        raise AttributeError(name) from exc


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_module.main())
