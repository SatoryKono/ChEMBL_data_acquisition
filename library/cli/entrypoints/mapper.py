"""Console script entry point for :mod:`library.utils.cli_tools.mapper_main`."""

from __future__ import annotations

from collections.abc import Sequence

from ._legacy import execute


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to the legacy ``mapper`` implementation."""

    return execute("mapper", argv)


__all__ = ["main"]

