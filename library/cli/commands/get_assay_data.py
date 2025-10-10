from __future__ import annotations

from collections.abc import Sequence

from library.cli.entrypoints import assay as _assay_entrypoint


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the assay pipeline CLI."""

    return _assay_entrypoint.main(argv)


__all__ = ["main"]
