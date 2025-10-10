"""Compatibility wrapper exposing :mod:`library.cli.commands.get_activities`."""

from __future__ import annotations

import sys
from collections.abc import Sequence
from importlib import import_module
from types import ModuleType

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a package module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

_MODULE: ModuleType = import_module("library.cli.commands.get_activities")
_ORIGINAL_MAIN = _MODULE.main
__all__ = tuple(getattr(_MODULE, "__all__", ()))


def __getattr__(name: str) -> object:
    """Delegate attribute access to :mod:`library.cli.commands.get_activities`."""

    try:
        return getattr(_MODULE, name)
    except AttributeError as exc:  # pragma: no cover - propagate missing attrs
        raise AttributeError(name) from exc


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the ``get-activities`` helper via the command module."""

    if argv is None:
        argv = list(sys.argv[1:])
    return _ORIGINAL_MAIN(argv)


sys.modules.setdefault("scripts.get_activities", sys.modules[__name__])


if __name__ == "__main__":
    raise SystemExit(main())
