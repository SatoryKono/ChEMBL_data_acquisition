"""Console script entry point for :mod:`library.utils.cli_tools.table_quality_main`."""

from __future__ import annotations

from collections.abc import Sequence

from ._legacy import execute


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to the legacy ``table_quality`` implementation."""

    return execute("table_quality", argv)


__all__ = ["main"]

