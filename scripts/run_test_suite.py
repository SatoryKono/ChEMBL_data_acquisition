"""Compatibility wrapper delegating to :mod:`scripts.run_tests`."""

from __future__ import annotations

import argparse
import os
import sys
from collections.abc import Sequence
from pathlib import Path

from scripts import run_tests

_DEPRECATION_MESSAGE = (
    "scripts.run_test_suite is deprecated and will be removed in a future "
    "release. Use `python scripts/run_tests.py` instead."
)

_DEFAULT_JSON_NAME = "test_report.json"
_DEFAULT_MARKDOWN_NAME = "test_summary.md"


def _preprocess_argv(argv: Sequence[str] | None) -> list[str]:
    if not argv:
        return []
    processed: list[str] = []
    skip_dashdash = False
    for token in argv:
        if token == "--pytest-args":
            skip_dashdash = True
            processed.append(token)
            continue
        if skip_dashdash and token == "--":
            skip_dashdash = False
            continue
        processed.append(token)
    return processed


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Execute pytest with reporting helpers (deprecated wrapper)."
    )
    parser.add_argument("--suite", default="full", help="Ignored legacy option.")
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable DEBUG logging in the delegated runner.",
    )
    parser.add_argument(
        "--report-dir",
        default="reports",
        type=Path,
        help="Directory where JSON and Markdown summaries will be written.",
    )
    parser.add_argument(
        "--pytest-args",
        nargs=argparse.REMAINDER,
        help="Additional pytest arguments forwarded after '--'.",
    )
    return parser


def _normalise_extra_args(extra_args: Sequence[str] | None) -> list[str]:
    if not extra_args:
        return []
    if extra_args and extra_args[0] == "--":
        return list(extra_args[1:])
    return list(extra_args)


def _emit_warning(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr, flush=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(_preprocess_argv(argv))

    os.environ.setdefault("PYTHONHASHSEED", "0")

    report_dir = Path(args.report_dir)
    json_path = report_dir / _DEFAULT_JSON_NAME
    markdown_path = report_dir / _DEFAULT_MARKDOWN_NAME

    run_args: list[str] = ["--json", str(json_path), "--markdown", str(markdown_path)]
    if args.verbose:
        run_args.insert(0, "--verbose")

    extra_args = _normalise_extra_args(args.pytest_args)
    if extra_args:
        run_args.append("--")
        run_args.extend(extra_args)

    _emit_warning(_DEPRECATION_MESSAGE)
    if args.suite and args.suite != "full":
        _emit_warning(f"Ignoring legacy --suite value: {args.suite}")

    return int(run_tests.main(run_args))


if __name__ == "__main__":
    raise SystemExit(main())
