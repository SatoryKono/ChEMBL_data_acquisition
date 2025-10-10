"""Console script entry point for :mod:`scripts.check_determinism`."""

from __future__ import annotations

from collections.abc import Sequence

from ._legacy import execute


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to the legacy ``check_determinism`` implementation."""

    return execute("check_determinism", argv)


__all__ = ["main"]

