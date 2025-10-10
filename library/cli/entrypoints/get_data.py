"""Console script entry point for :mod:`scripts.get_data`."""

from __future__ import annotations

from collections.abc import Sequence

from ._legacy import execute


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to the legacy ``get_data`` implementation."""

    return execute("get_data", argv)


__all__ = ["main"]

