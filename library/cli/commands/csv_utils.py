from __future__ import annotations

from collections.abc import Sequence

from . import _run


def main(argv: Sequence[str] | None = None) -> int:
    """Run scripts.csv_utils.

    Parameters
    ----------
    argv:
        Optional sequence of command-line arguments.

    Returns
    -------
    int
        Exit code returned by the script.
    """

    return _run("csv_utils", argv)
