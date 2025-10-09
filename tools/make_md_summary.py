"""Render a Markdown summary from a structured pytest JSON report."""
from __future__ import annotations

# ruff: noqa: E402  # requires sys.path mutation before local imports

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:  # pragma: no cover - import side effect
    sys.path.insert(0, str(ROOT_DIR))

# ruff: noqa: E402  # allow importing project modules after sys.path adjustments
from scripts.run_tests import build_summary_markdown

DEFAULT_INPUT = Path("reports/test_report.json")
DEFAULT_OUTPUT = Path("reports/test_summary.md")
REQUIRED_SUMMARY_KEYS = {
    "total",
    "passed",
    "failed",
    "skipped",
    "xfailed",
    "xpassed",
    "error",
    "success_rate",
}


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render a Markdown summary from a structured pytest JSON report."
        )
    )
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        default=DEFAULT_INPUT,
        help=(
            "Path to the structured JSON report produced by pytest. "
            "Defaults to reports/test_report.json."
        ),
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=(
            "Destination path for the Markdown summary. "
            "Defaults to reports/test_summary.md."
        ),
    )
    return parser.parse_args(argv)


def _load_report(path: Path) -> dict[str, Any]:
    try:
        text = path.read_text(encoding="utf-8")
    except FileNotFoundError as exc:  # pragma: no cover - defensive guard
        raise RuntimeError(f"Input report {path} does not exist") from exc
    except OSError as exc:  # pragma: no cover - filesystem issues
        raise RuntimeError(f"Unable to read input report {path}: {exc}") from exc

    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Input report {path} is not valid JSON: {exc}") from exc


def _write_summary(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _validate_report(report: dict[str, Any]) -> None:
    if not isinstance(report, dict):
        raise ValueError("Report payload must be a JSON object")

    summary = report.get("summary")
    if not isinstance(summary, dict):
        raise ValueError("Report missing 'summary' section")

    missing = REQUIRED_SUMMARY_KEYS - summary.keys()
    if missing:
        missing_list = ", ".join(sorted(missing))
        raise ValueError(f"Summary missing required keys: {missing_list}")

    for key in REQUIRED_SUMMARY_KEYS - {"success_rate"}:
        value = summary.get(key)
        if not isinstance(value, int):
            raise ValueError(
                f"Summary field '{key}' must be an integer, got {type(value).__name__}"
            )

    success_rate = summary.get("success_rate")
    if not isinstance(success_rate, (int, float)):
        raise ValueError("Summary field 'success_rate' must be numeric")

    tests = report.get("tests")
    if not isinstance(tests, list):
        raise ValueError("Report 'tests' section must be a list")
    for index, entry in enumerate(tests):
        if not isinstance(entry, dict):
            raise ValueError(
                f"Test entry at index {index} must be an object, "
                f"got {type(entry).__name__}"
            )

    meta = report.get("meta")
    if meta is not None and not isinstance(meta, dict):
        raise ValueError("Report 'meta' section must be an object when present")


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)

    try:
        report = _load_report(args.input)
        _validate_report(report)
        summary_markdown = build_summary_markdown(report)
        _write_summary(args.output, summary_markdown)
    except ValueError as exc:
        print(f"Invalid report: {exc}", file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print(exc, file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entrypoint
    raise SystemExit(main())
