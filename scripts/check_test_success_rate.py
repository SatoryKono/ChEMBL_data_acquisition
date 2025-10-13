"""Utility to validate pytest success rate from structured JSON reports."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

ROOT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_REPORT_PATH = ROOT_DIR / "reports" / "test_report.json"
DEFAULT_THRESHOLD = 95.0


class ReportValidationError(RuntimeError):
    """Raised when the report payload is missing required information."""


def _resolve_report_path(raw_path: str | None) -> Path:
    if not raw_path:
        return DEFAULT_REPORT_PATH
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = ROOT_DIR / candidate
    return candidate


def _load_report(path: Path) -> dict[str, Any]:
    try:
        raw_text = path.read_text(encoding="utf-8")
    except FileNotFoundError as exc:  # pragma: no cover - exercised in CI
        raise ReportValidationError(
            f"Report file '{path}' does not exist; ensure tests generated it"
        ) from exc
    try:
        data = json.loads(raw_text)
    except json.JSONDecodeError as exc:
        raise ReportValidationError(
            f"Report file '{path}' contains invalid JSON: {exc}"
        ) from exc
    if not isinstance(data, dict):
        raise ReportValidationError("Report payload must be a JSON object")
    return data


def _extract_success_rate(summary: dict[str, Any]) -> float:
    raw_rate = summary.get("success_rate")
    if isinstance(raw_rate, int | float):
        rate = float(raw_rate)
        if rate > 1.0:
            rate = rate / 100.0
        if 0.0 <= rate <= 1.0:
            return rate
    total_raw = summary.get("total")
    passed_raw = summary.get("passed")
    skipped_raw = summary.get("skipped", 0)
    xfailed_raw = summary.get("xfailed", 0)

    field_names = {
        "total": total_raw,
        "passed": passed_raw,
        "skipped": skipped_raw,
        "xfailed": xfailed_raw,
    }
    validated: dict[str, int] = {}
    for name, value in field_names.items():
        if not isinstance(value, int):
            raise ReportValidationError(
                f"Summary field '{name}' must be an integer to compute success rate"
            )
        if value < 0:
            raise ReportValidationError(f"Summary field '{name}' must be non-negative")
        validated[name] = value

    total = validated["total"]
    passed = validated["passed"]
    skipped = validated["skipped"]
    xfailed = validated["xfailed"]

    denominator = max(1, total - skipped)
    numerator = passed + xfailed
    return max(0.0, min(1.0, numerator / denominator))


def _parse_threshold(value: float) -> float:
    if value < 0 or value > 100:
        raise argparse.ArgumentTypeError(
            "Threshold must be a percentage between 0 and 100"
        )
    return value


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate pytest success rate from a structured JSON report",
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Path to reports/test_report.json (relative to repository root)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help="Minimum required success rate in percent (default: 95.0)",
    )
    parser.add_argument(
        "--require-total",
        type=int,
        default=None,
        help="Optional minimum number of collected tests to consider the report valid",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_argument_parser()
    args = parser.parse_args(argv)

    report_path = _resolve_report_path(args.report)
    try:
        threshold_pct = _parse_threshold(float(args.threshold))
    except argparse.ArgumentTypeError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    threshold_rate = threshold_pct / 100.0

    try:
        payload = _load_report(report_path)
    except ReportValidationError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    summary = payload.get("summary")
    if not isinstance(summary, dict):
        print("Report summary section is missing or malformed", file=sys.stderr)
        return 2

    if args.require_total is not None:
        required_total = int(args.require_total)
        actual_total = summary.get("total")
        if not isinstance(actual_total, int) or actual_total < required_total:
            print(
                "Report total tests ({actual}) is below the required minimum ({required})".format(
                    actual=(
                        actual_total if isinstance(actual_total, int) else "<unknown>"
                    ),
                    required=required_total,
                ),
                file=sys.stderr,
            )
            return 2

    try:
        success_rate = _extract_success_rate(summary)
    except ReportValidationError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    success_pct = success_rate * 100.0

    total = summary.get("total", "<unknown>")
    passed = summary.get("passed", "<unknown>")

    if success_rate + 1e-9 < threshold_rate:
        print(
            f"Test success rate {success_pct:.2f}% is below the required threshold {threshold_pct:.2f}% (passed {passed} / total {total})",
            file=sys.stderr,
        )
        return 1

    print(
        f"Test success rate {success_pct:.2f}% meets the required threshold {threshold_pct:.2f}% (passed {passed} / total {total})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
