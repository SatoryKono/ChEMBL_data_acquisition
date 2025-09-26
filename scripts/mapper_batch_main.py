"""Proxy module for :mod:`library.utils.cli_tools.mapper_batch_main`."""

from __future__ import annotations

from collections.abc import Sequence

from library.utils.cli_tools import mapper_batch_main as _impl


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to :func:`library.utils.cli_tools.mapper_batch_main.main`."""

    return _impl.main(argv)


if __name__ == "__main__":  # pragma: no cover - convenience entry point
    raise SystemExit(main())
