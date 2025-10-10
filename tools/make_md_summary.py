"""Render a Markdown summary from a structured pytest JSON report."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any


# Ensure the repository root is importable when the script is executed directly
# via ``python tools/make_md_summary.py`` from a fresh checkout.
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


from library.reporting.test_summary import (  # noqa: E402  - depends on path tweak above
    build_summary_markdown,
    validate_summary_report,
)

DEFAULT_INPUT = Path("reports/test_report.json")
DEFAULT_OUTPUT = Path("reports/test_summary.md")


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Render a Markdown summary from a structured pytest JSON report.")
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


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)

    try:
        report = _load_report(args.input)
        validate_summary_report(report)
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
