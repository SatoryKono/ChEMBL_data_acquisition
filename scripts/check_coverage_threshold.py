"""Validate that pytest coverage meets the configured threshold."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from xml.etree import ElementTree

ROOT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_COVERAGE_XML = ROOT_DIR / "reports" / "coverage" / "coverage.xml"
DEFAULT_THRESHOLD = 80.0


class CoverageValidationError(RuntimeError):
    """Raised when the coverage report is missing or malformed."""


def _resolve_coverage_path(raw_path: str | None) -> Path:
    if not raw_path:
        return DEFAULT_COVERAGE_XML
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = ROOT_DIR / candidate
    return candidate


def _load_line_rate(path: Path) -> float:
    try:
        xml_text = path.read_text(encoding="utf-8")
    except FileNotFoundError as exc:  # pragma: no cover - used in CI when missing
        raise CoverageValidationError(
            f"Coverage report '{path}' does not exist; ensure pytest generated it"
        ) from exc

    try:
        root = ElementTree.fromstring(xml_text)
    except ElementTree.ParseError as exc:
        raise CoverageValidationError(
            f"Coverage report '{path}' contains invalid XML: {exc}"
        ) from exc

    line_rate_raw = root.attrib.get("line-rate")
    if line_rate_raw is None:
        raise CoverageValidationError(
            "Coverage XML is missing the required 'line-rate' attribute"
        )

    try:
        line_rate = float(line_rate_raw)
    except ValueError as exc:
        raise CoverageValidationError(
            "Coverage XML attribute 'line-rate' must be numeric"
        ) from exc

    if line_rate > 1.0:
        return max(0.0, min(100.0, line_rate))
    return max(0.0, min(100.0, line_rate * 100.0))


def _parse_threshold(value: float) -> float:
    if not 0.0 <= value <= 100.0:
        raise argparse.ArgumentTypeError(
            "Threshold must be a percentage between 0 and 100"
        )
    return value


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate that line coverage reported by coverage.py meets the threshold",
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Path to coverage.xml (relative to repository root)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help="Minimum required line coverage percentage (default: 80.0)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    report_path = _resolve_coverage_path(args.report)
    try:
        threshold_pct = _parse_threshold(float(args.threshold))
    except argparse.ArgumentTypeError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    try:
        coverage_pct = _load_line_rate(report_path)
    except CoverageValidationError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    if coverage_pct + 1e-9 < threshold_pct:
        print(
            f"Line coverage {coverage_pct:.2f}% is below the required threshold {threshold_pct:.2f}%",
            file=sys.stderr,
        )
        return 1

    print(
        f"Line coverage {coverage_pct:.2f}% meets the required threshold {threshold_pct:.2f}%",
    )
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
