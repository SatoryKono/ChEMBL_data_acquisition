"""Console script entry point for :mod:`library.utils.cli_tools.get_input_initialisation`."""

from __future__ import annotations

from collections.abc import Sequence

from ._legacy import execute


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to the legacy ``get_input_initialisation`` implementation."""

    return execute("get_input_initialisation", argv)


__all__ = ["main"]

