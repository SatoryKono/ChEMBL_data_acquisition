"""Compatibility wrapper exposing the activity pipeline CLI entry point."""

from __future__ import annotations

import sys
from importlib import import_module
from types import ModuleType

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli


def _load_activity_module() -> ModuleType:
    module = import_module("library.cli.commands.get_activity_data")
    return module


_activity = _load_activity_module()
__all__ = tuple(getattr(_activity, "__all__", ()))


def __getattr__(name: str) -> object:
    """Proxy attribute access to :mod:`library.cli.commands.get_activity_data`."""

    try:
        return getattr(_activity, name)
    except AttributeError as exc:  # pragma: no cover - passthrough for missing attrs
        raise AttributeError(name) from exc


sys.modules[__name__] = _activity
sys.modules.setdefault("scripts.get_activity_data", _activity)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_activity.main())
