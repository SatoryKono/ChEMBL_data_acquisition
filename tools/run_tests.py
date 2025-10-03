"""Execute pytest and emit JSON and Markdown reports."""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pytest

from scripts.run_test_suite import SUCCESS_RATE_THRESHOLD


REPO_NAME = "SatoryKono/ChEMBL_data_acquisition"


logger = logging.getLogger(__name__)

@dataclass
class TestRecord:
    nodeid: str
    status: str
    duration_ms: float
    stdout: str
    stderr: str
    log: list[dict[str, str]] = field(default_factory=list)
    error: str | None = None


class ReportCollector:
    def __init__(self) -> None:
        self.tests: dict[str, TestRecord] = {}
        self.summary = {
            "total": 0,
            "passed": 0,
            "failed": 0,
            "skipped": 0,
            "xfailed": 0,
            "xpassed": 0,
            "error": 0,
        }
        self.start_time: float = 0.0
        self.end_time: float = 0.0

    # pytest hooks -----------------------------------------------------
    def pytest_sessionstart(self, session: pytest.Session) -> None:  # noqa: D401
        self.start_time = time.perf_counter()

    def pytest_sessionfinish(self, session: pytest.Session) -> None:  # noqa: D401
        self.end_time = time.perf_counter()

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:  # noqa: D401
        if report.when not in {"call", "setup", "teardown"}:
            return
        key = report.nodeid
        record = self.tests.get(key)
        if record is None:
            record = TestRecord(
                nodeid=key,
                status="",
                duration_ms=0.0,
                stdout="",
                stderr="",
            )
            self.tests[key] = record

        if report.when == "call" or (report.failed and record.status == ""):
            status = self._determine_status(report)
            record.status = status
            record.duration_ms = max(report.duration * 1000.0, 0.0)
            record.stdout = getattr(report, "capstdout", "") or ""
            record.stderr = getattr(report, "capstderr", "") or ""
            record.log = [
                {"section": name, "content": content}
                for name, content in getattr(report, "sections", [])
            ]
            if report.failed:
                record.error = self._format_longrepr(report.longrepr)
        elif record.status == "":
            # Handle errors during setup/teardown
            status = self._determine_status(report)
            record.status = status
            record.duration_ms = max(record.duration_ms, report.duration * 1000.0)
            if report.failed:
                record.error = self._format_longrepr(report.longrepr)

    # helper methods ---------------------------------------------------
    @staticmethod
    def _determine_status(report: pytest.TestReport) -> str:
        if report.skipped:
            if getattr(report, "wasxfail", False):
                return "xfailed"
            return "skipped"
        if report.failed:
            return "failed" if report.when == "call" else "error"
        if getattr(report, "wasxfail", False):
            return "xpassed"
        return "passed"

    @staticmethod
    def _format_longrepr(longrepr: Any) -> str:
        if longrepr is None:
            return ""
        if hasattr(longrepr, "reprcrash"):
            return str(longrepr)
        return str(longrepr)

    # synthesis --------------------------------------------------------
    def build_results(self) -> list[TestRecord]:
        results = list(self.tests.values())
        results.sort(key=lambda record: record.nodeid)
        self.summary["total"] = len(results)
        for record in results:
            status = record.status or "error"
            if status in self.summary:
                self.summary[status] += 1
        return results

    def duration(self) -> float:
        end = self.end_time or time.perf_counter()
        start = self.start_time
        return max(end - start, 0.0)


def _git_output(args: list[str]) -> str:
    try:
        completed = subprocess.run(
            ["git", *args],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError:
        return ""
    return completed.stdout.strip()


def _build_json(
    collector: ReportCollector,
    *,
    report_path: Path,
) -> dict[str, Any]:
    results = collector.build_results()
    summary = collector.summary
    duration_sec = collector.duration()
    if summary["total"]:
        success_rate = round(100.0 * summary["passed"] / summary["total"], 2)
    else:
        success_rate = 100.0

    repo = _git_output(["config", "--get", "remote.origin.url"]) or REPO_NAME
    payload = {
        "meta": {
            "repo": repo,
            "commit": _git_output(["rev-parse", "HEAD"]),
            "branch": _git_output(["rev-parse", "--abbrev-ref", "HEAD"]),
            "ts_utc": datetime.now(timezone.utc).isoformat(),
            "duration_sec": round(duration_sec, 3),
            "python": sys.version.split()[0],
            "pytest": pytest.__version__,
        },
        "summary": {
            **summary,
            "success_rate": success_rate,
        },
        "tests": [
            {
                "nodeid": record.nodeid,
                "status": record.status,
                "duration_ms": round(record.duration_ms, 3),
                "stdout": record.stdout,
                "stderr": record.stderr,
                "log": record.log,
                "error": record.error,
            }
            for record in results
        ],
    }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def _write_markdown(summary_path: Path, data: dict[str, Any]) -> None:
    summary = data["summary"]
    meta = data["meta"]
    lines = [
        "# Test Summary",
        "",
        f"- Repo: `{meta['repo']}`",
        f"- Commit: {meta['commit']}",
        f"- Branch: {meta['branch']}",
        f"- Timestamp (UTC): {meta['ts_utc']}",
        f"- Duration: {meta['duration_sec']} s",
        f"- Success rate: {summary['success_rate']}%",
        "",
        "| total | passed | failed | skipped | xfailed | xpassed | error |",
        "|------:|-------:|-------:|--------:|--------:|--------:|------:|",
        "| {total:5d} | {passed:5d} | {failed:5d} | {skipped:6d} | {xfailed:6d} | {xpassed:6d} | {error:5d} |".format(
            total=summary["total"],
            passed=summary["passed"],
            failed=summary["failed"],
            skipped=summary["skipped"],
            xfailed=summary["xfailed"],
            xpassed=summary["xpassed"],
            error=summary["error"],
        ),
        "",
    ]

    failed_entries = [
        (test["nodeid"], test["error"] or test["stderr"] or "")
        for test in data["tests"]
        if test["status"] in {"failed", "error"}
    ]
    if failed_entries:
        lines.append("## Failed / Error details")
        for nodeid, message in failed_entries:
            short = message.splitlines()[0] if message else "(no message)"
            lines.append(f"- `{nodeid}`: `{short}`")
    else:
        lines.append("## Failed / Error details")
        lines.append("- none")

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run pytest and emit reports")
    parser.add_argument(
        "--pytest-args",
        nargs=argparse.REMAINDER,
        help="Additional arguments forwarded to pytest",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=Path("reports/test_report.json"),
        help="Path to JSON report output",
    )
    parser.add_argument(
        "--markdown",
        type=Path,
        default=Path("reports/test_summary.md"),
        help="Path to Markdown summary",
    )
    args = parser.parse_args(argv)

    collector = ReportCollector()
    pytest_args = ["-q"]
    if args.pytest_args:
        pytest_args.extend(args.pytest_args)

    pytest_exit_code = pytest.main(pytest_args, plugins=[collector])
    exit_code = int(pytest_exit_code)

    data = _build_json(collector, report_path=args.json)
    _write_markdown(args.markdown, data)

    success_rate_percent = data["summary"]["success_rate"]
    success_rate_ratio = success_rate_percent / 100.0
    if success_rate_ratio < SUCCESS_RATE_THRESHOLD:
        logger.error(
            "Success rate %.2f%% is below the required threshold of %.2f%%",
            success_rate_percent,
            SUCCESS_RATE_THRESHOLD * 100,
        )
        if exit_code == 0:
            exit_code = 1

    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())

