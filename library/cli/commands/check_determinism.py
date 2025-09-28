from __future__ import annotations

from collections.abc import Sequence

from . import _run


def main(argv: Sequence[str] | None = None) -> int:
    """Run :mod:`library.utils.cli_tools.check_determinism`.

    Parameters
    ----------
    argv:
        Optional sequence of command-line arguments.

    Returns
    -------
    int
        Exit code returned by the script.
    """

    return _run("check_determinism", argv)
