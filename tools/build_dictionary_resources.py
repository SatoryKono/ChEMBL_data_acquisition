"""Utility entrypoint to rebuild bundled dictionary resources.

The helper delegates to existing CLI modules that know how to refresh
individual caches.  It does not aim to replicate their command line
interfaces; instead it provides a documented place referenced by the
manifest so that engineers can discover the appropriate regeneration
commands.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

_TARGET_SCRIPT = PROJECT_ROOT / "scripts" / "get_target_data.py"
_TESTITEM_SCRIPT = PROJECT_ROOT / "scripts" / "get_testitem_data.py"


def _run(command: list[str]) -> None:
    """Run *command* printing a concise log message on failure."""

    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:  # pragma: no cover - helper
        raise SystemExit(exc.returncode) from exc


def _build_target_resources() -> None:
    """Recompute target dictionaries shipped with the repository."""

    _run([
        "python",
        str(_TARGET_SCRIPT),
        "--mode",
        "all",
        "--dry-run",
    ])


def _build_testitem_resources() -> None:
    """Recompute test item lookup tables shipped with the repository."""

    _run([
        "python",
        str(_TESTITEM_SCRIPT),
        "--dry-run",
    ])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--targets",
        action="store_true",
        help="Regenerate target dictionaries using scripts/get_target_data.py",
    )
    parser.add_argument(
        "--testitems",
        action="store_true",
        help="Regenerate test item dictionaries using scripts/get_testitem_data.py",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Regenerate both target and test item dictionaries",
    )
    args = parser.parse_args(argv)

    if not any((args.targets, args.testitems, args.all)):
        parser.print_help()
        return 0

    if args.targets or args.all:
        _build_target_resources()
    if args.testitems or args.all:
        _build_testitem_resources()
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI wrapper
    raise SystemExit(main())
