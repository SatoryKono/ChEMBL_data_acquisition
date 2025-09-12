from __future__ import annotations

from collections.abc import Sequence

from . import _run


def main(argv: Sequence[str] | None = None) -> int:
    """Run scripts.chunk_io.

    Parameters
    ----------
    argv:
        Optional sequence of command-line arguments.

    Returns
    -------
    int
        Exit code returned by the script.
    """

    return _run("chunk_io", argv)
