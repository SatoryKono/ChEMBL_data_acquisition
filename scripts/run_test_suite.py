"""Run the full pytest suite and export structured reports."""

from __future__ import annotations

import argparse
import json
import logging
import os
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence

import pytest


from library.common.logging_setup import LoggerConfig, configure_logger
from library.cli.logging import setup_cli_logging


SUCCESS_RATE_THRESHOLD = 0.95


logger = logging.getLogger(__name__)

ROOT_DIR = Path(__file__).resolve().parents[1]
TEST_DIRECTORIES = (
    ROOT_DIR / "tests" / "unit",
    ROOT_DIR / "tests" / "integration",
    ROOT_DIR / "tests" / "postprocessing",
    ROOT_DIR / "tests" / "e2e",
)


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


def summarize_results(results: Sequence[TestResult]) -> dict[str, Any]:
    total = len(results)
    passed = sum(1 for item in results if item.status == "passed")
    failed = sum(1 for item in results if item.status == "failed")
    skipped = sum(1 for item in results if item.status == "skipped")
    duration = sum(item.duration for item in results)
    success_rate = 1.0 if total == 0 else passed / total

    return {
        "total": total,
        "passed": passed,
        "failed": failed,
        "skipped": skipped,
        "duration": duration,
        "success_rate": success_rate,
    }


def _write_json(report_path: Path, results: Sequence[TestResult], exit_code: int, summary: dict[str, Any]) -> None:
    payload = {
        "generated_at": datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "exit_code": exit_code,
        "summary": summary,
        "tests": [asdict(result) for result in results],
    }
    report_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _write_summary(
    summary_path: Path,
    results: Sequence[TestResult],
    exit_code: int,
    summary: dict[str, Any],
) -> None:
    total = summary["total"]
    passed = summary["passed"]
    failed = summary["failed"]
    skipped = summary["skipped"]
    duration = summary["duration"]
    success_rate = summary["success_rate"] * 100

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
        f"* Success rate: {success_rate:.2f}%",
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
        "--verbose",
        action="store_true",
        help="Enable debug logging for the CLI and pytest log capture.",
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
        help="Additional arguments forwarded to pytest after '--'.",
    )
    return parser


def _has_explicit_targets(extra_args: Sequence[str] | None) -> bool:
    if not extra_args:
        return False
    for token in extra_args:
        if token == "--":
            continue
        if token.startswith("-"):
            continue
        return True
    return False


def _normalise_extra_args(extra_args: Sequence[str] | None) -> list[str]:
    if not extra_args:
        return []
    if extra_args and extra_args[0] == "--":
        return list(extra_args[1:])
    return list(extra_args)


def _default_test_targets() -> list[str]:
    return [str(path) for path in TEST_DIRECTORIES if path.exists()]


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    report_dir: Path = args.report_dir
    report_dir.mkdir(parents=True, exist_ok=True)

    if "PYTHONHASHSEED" not in os.environ:
        os.environ["PYTHONHASHSEED"] = "0"

    level = "DEBUG" if args.verbose else "INFO"
    base_logger_cfg = LoggerConfig(level=level, logger_name="run_test_suite")

    with setup_cli_logging("run_test_suite", base_logger_cfg) as logging_ctx:
        configure_logger(logging_ctx.log_cfg)

        log_path = logging_ctx.log_path
        extra_args = _normalise_extra_args(args.pytest_args)

        pytest_args = [
            f"--log-file={log_path}",
            f"--log-file-level={logging_ctx.log_cfg.level}",
            "--maxfail=0",
        ]

        if not _has_explicit_targets(extra_args):
            pytest_args.extend(_default_test_targets())

        pytest_args.extend(extra_args)

        plugin = JsonReportPlugin(log_path)
        pytest_exit_code = pytest.main(pytest_args, plugins=[plugin])
        exit_code = int(pytest_exit_code)

        results = plugin.results
        summary = summarize_results(results)

        if summary["success_rate"] < SUCCESS_RATE_THRESHOLD:
            logger.error(
                "Success rate %.2f%% is below the required threshold of %.2f%%",
                summary["success_rate"] * 100,
                SUCCESS_RATE_THRESHOLD * 100,
            )
            if exit_code == 0:
                exit_code = 1

        json_path = report_dir / "test_report.json"
        summary_path = report_dir / "test_summary.md"
        _write_json(json_path, results, exit_code, summary)
        _write_summary(summary_path, results, exit_code, summary)

    return int(exit_code)


if __name__ == "__main__":
    raise SystemExit(main())
