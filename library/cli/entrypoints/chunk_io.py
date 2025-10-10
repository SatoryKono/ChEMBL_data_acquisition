"""Console script entry point for :mod:`library.utils.cli_tools.chunk_io_main`."""

from __future__ import annotations

from collections.abc import Sequence

from ._legacy import execute


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to the legacy ``chunk_io`` implementation."""

    return execute("chunk_io", argv)


__all__ = ["main"]

