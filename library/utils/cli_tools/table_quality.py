"""Proxy module for :mod:`library.utils.cli_tools.table_quality_main`."""

from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path

if __package__ in {None, ""}:
    project_root = Path(__file__).resolve().parents[3]
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

from library.utils import bootstrap
from library.utils.cli_tools import table_quality_main as _impl

bootstrap.ensure_project_root()


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to :func:`library.utils.cli_tools.table_quality_main.main`."""

    return _impl.main(argv)


if __name__ == "__main__":  # pragma: no cover - convenience entry point
    raise SystemExit(main())
