"""Helpers for delegating CLI entry points to legacy implementations."""

from __future__ import annotations

from collections.abc import Sequence

from library.cli.commands import _run as _run_legacy_command


def execute(module: str, argv: Sequence[str] | None = None) -> int:
    """Execute a legacy CLI module and return its exit code.

    Parameters
    ----------
    module:
        Name of the legacy module to execute.  The value is forwarded to
        :func:`library.cli.commands._run`, which resolves the module either in
        ``library.utils.cli_tools`` or ``scripts`` namespaces.
    argv:
        Optional command line arguments passed to the legacy module.
    """

    return _run_legacy_command(module, argv)


__all__ = ["execute"]

