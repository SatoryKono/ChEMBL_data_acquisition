"""Run the full pytest suite and export structured reports."""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Any

import pytest


@dataclass
class TestResult:
    """Serializable representation of a single pytest outcome."""

    name: str
    status: str
    duration: float
    message: str
    log_path: str | None


class JsonReportPlugin:
    """Collect pytest results and persist optional per-suite logs."""

    def __init__(self, log_path: Path) -> None:
        self._log_path = log_path
        self._results: dict[str, TestResult] = {}

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:  # type: ignore[override]
        entry = self._results.get(report.nodeid)
        if entry is None:
            entry = TestResult(
                name=report.nodeid,
                status="passed",
                duration=0.0,
                message="",
                log_path=str(self._log_path) if self._log_path else None,
            )
            self._results[report.nodeid] = entry

        entry.duration += report.duration

        if report.outcome == "failed":
            entry.status = "failed"
            entry.message = _format_longrepr(report.longrepr)
        elif report.outcome == "skipped":
            entry.status = "skipped"
            entry.message = _format_longrepr(report.longrepr)

    @property
    def results(self) -> list[TestResult]:
        return list(self._results.values())


def _format_longrepr(longrepr: Any) -> str:
    if isinstance(longrepr, tuple) and len(longrepr) == 3:
        return str(longrepr[2])
    if hasattr(longrepr, "longreprtext"):
        return str(longrepr.longreprtext)  # pragma: no cover - compatibility guard
    return str(longrepr)


def _write_json(report_path: Path, plugin: JsonReportPlugin, exit_code: int) -> None:
    payload = {
        "generated_at": datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "exit_code": exit_code,
        "tests": [asdict(result) for result in plugin.results],
    }
    report_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _write_summary(summary_path: Path, plugin: JsonReportPlugin, exit_code: int) -> None:
    results = plugin.results
    total = len(results)
    passed = sum(1 for item in results if item.status == "passed")
    failed = sum(1 for item in results if item.status == "failed")
    skipped = sum(1 for item in results if item.status == "skipped")
    duration = sum(item.duration for item in results)

    lines = [
        "# Test suite summary",
        "",
        f"* Generated at: {datetime.utcnow().isoformat(timespec='seconds')}Z",
        f"* Exit code: {exit_code}",
        f"* Tests collected: {total}",
        f"* Passed: {passed}",
        f"* Failed: {failed}",
        f"* Skipped: {skipped}",
        f"* Cumulative duration: {duration:.2f}s",
        "",
        "## Failing tests" if failed else "## Failing tests",
    ]

    if failed:
        for item in results:
            if item.status != "failed":
                continue
            lines.extend(
                [
                    f"- `{item.name}` ({item.duration:.2f}s)",
                    f"  - Log: `{item.log_path}`" if item.log_path else "  - Log: not captured",
                    "  - Message:",
                    "    ```",
                    *[f"    {line}" for line in item.message.splitlines() or ["<no message>"]],
                    "    ```",
                ]
            )
    else:
        lines.append("- None")

    skipped_section_header = "## Skipped tests"
    lines.append("")
    lines.append(skipped_section_header)
    if skipped:
        for item in results:
            if item.status != "skipped":
                continue
            lines.extend(
                [
                    f"- `{item.name}` ({item.duration:.2f}s)",
                    f"  - Reason: {item.message or 'not provided'}",
                ]
            )
    else:
        lines.append("- None")

    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Execute pytest with reporting helpers.")
    parser.add_argument("--suite", default="full", help="Label for the executed suite (used in log naming).")
    parser.add_argument(
        "--report-dir",
        default="reports",
        type=Path,
        help="Directory where JSON and Markdown summaries will be written.",
    )
    parser.add_argument(
        "--pytest-args",
        nargs=argparse.REMAINDER,
        help="Additional arguments forwarded to pytest after '--'.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    report_dir: Path = args.report_dir
    report_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = report_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    log_path = logs_dir / f"{args.suite}.log"

    if "PYTHONHASHSEED" not in os.environ:
        os.environ["PYTHONHASHSEED"] = "0"

    pytest_args = ["tests", f"--log-file={log_path}", "--log-file-level=DEBUG", "--maxfail=0"]
    if args.pytest_args:
        pytest_args.extend(args.pytest_args)

    plugin = JsonReportPlugin(log_path)
    exit_code = pytest.main(pytest_args, plugins=[plugin])

    json_path = report_dir / "test_report.json"
    summary_path = report_dir / "test_summary.md"
    _write_json(json_path, plugin, exit_code)
    _write_summary(summary_path, plugin, exit_code)

    return int(exit_code)


if __name__ == "__main__":
    raise SystemExit(main())
