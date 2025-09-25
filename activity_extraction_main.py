"""Backward-compatible CLI wrapper for activity extraction."""

from __future__ import annotations

from collections.abc import Sequence

from scripts import get_activity_data as _activity_cli

build_parser = _activity_cli.build_parser


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to :mod:`scripts.get_activity_data`."""

    return _activity_cli.main(argv)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
