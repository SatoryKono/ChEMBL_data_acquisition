from __future__ import annotations

from collections.abc import Sequence

from library.cli.entrypoints import document as _document_entrypoint


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the document pipeline CLI."""

    return _document_entrypoint.main(argv)


__all__ = ["main"]
