"""Shared helpers for CLI command wrappers."""

from __future__ import annotations

from collections.abc import Sequence
from importlib import import_module
from typing import cast


def _run(module: str, argv: Sequence[str] | None = None) -> int:
    """Execute ``module.main`` from the :mod:`scripts` package.

    Parameters
    ----------
    module:
        Name of the module inside :mod:`scripts`.
    argv:
        Optional sequence of command-line arguments.

    Returns
    -------
    int
        Exit code returned by the script's ``main`` function.
    """

    main_func = import_module(f"scripts.{module}").main
    return cast(int, main_func(argv))


__all__ = ["_run"]
