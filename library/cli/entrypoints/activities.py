"""Console script entry point for :mod:`scripts.get_activities`."""

from __future__ import annotations

from collections.abc import Sequence

from ._legacy import execute


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to the legacy ``get_activities`` implementation."""

    return execute("get_activities", argv)


__all__ = ["main"]

