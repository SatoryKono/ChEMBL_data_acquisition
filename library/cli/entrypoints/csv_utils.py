"""Console script entry point for :mod:`library.utils.cli_tools.csv_utils`."""

from __future__ import annotations

from collections.abc import Sequence

from ._legacy import execute


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to the legacy ``csv_utils`` implementation."""

    return execute("csv_utils", argv)


__all__ = ["main"]

