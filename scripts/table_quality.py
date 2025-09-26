"""Proxy module for :mod:`library.utils.cli_tools.table_quality_main`."""

from __future__ import annotations

from collections.abc import Sequence

from library.utils.cli_tools import table_quality_main as _impl


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to :func:`library.utils.cli_tools.table_quality_main.main`."""

    return _impl.main(argv)


if __name__ == "__main__":  # pragma: no cover - convenience entry point
    raise SystemExit(main())
