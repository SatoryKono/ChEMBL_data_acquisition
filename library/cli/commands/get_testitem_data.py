from __future__ import annotations

from collections.abc import Sequence

from library.cli.entrypoints import testitem as _testitem_entrypoint


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the test item pipeline CLI."""

    return _testitem_entrypoint.main(argv)


__all__ = ["main"]
