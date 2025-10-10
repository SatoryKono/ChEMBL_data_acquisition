"""Console script entry point for :mod:`library.utils.cli_tools.get_document_type`."""

from __future__ import annotations

from collections.abc import Sequence

from ._legacy import execute


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to the legacy ``get_document_type`` implementation."""

    return execute("get_document_type", argv)


__all__ = ["main"]

