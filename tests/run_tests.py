"""Execute pytest and emit JSON and Markdown reports."""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# ruff: noqa: E402  # the test harness adjusts sys.path to import project modules
import pytest

from library.cli.logging import setup_cli_logging
from library.common.logging_setup import LoggerConfig, configure_logger

REPO_NAME = "SatoryKono/ChEMBL_data_acquisition"
TEST_ROOT = Path(__file__).resolve().parent
TEST_DIRECTORIES = (
    TEST_ROOT / "unit",
    TEST_ROOT / "integration",
    TEST_ROOT / "postprocessing",
    TEST_ROOT / "e2e",
)

SUCCESS_RATE_THRESHOLD = 0.95


logger = logging.getLogger(__name__)

DEPRECATION_MESSAGE = (
    "tests/run_tests.py is deprecated and will be removed in a future release. "
    "Use scripts/run_tests.py instead."
)


@dataclass
class TestRecord:
    nodeid: str
    status: str
    duration_ms: float
    stdout: str
    stderr: str
    log: list[dict[str, str]] = field(default_factory=list)
    error: str | None = None
    artifacts: dict[str, object] = field(default_factory=dict)


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

        for key, value in getattr(report, "user_properties", ()):
            record.artifacts[key] = self._normalise_property_value(value)

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

    @staticmethod
    def _normalise_property_value(value: object) -> object:
        if isinstance(value, Path):
            return str(value)
        return value

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
            capture_output=True,
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

    executed_total = summary["total"] - summary["skipped"]
    executed_total = max(executed_total, 0)
    successes = summary["passed"] + summary["xfailed"]

    if executed_total:
        success_rate = successes / executed_total
    else:
        success_rate = 1.0

    success_rate = max(0.0, min(success_rate, 1.0))

    repo = _git_output(["config", "--get", "remote.origin.url"]) or REPO_NAME
    payload = {
        "meta": {
            "repo": repo,
            "commit": _git_output(["rev-parse", "HEAD"]),
            "branch": _git_output(["rev-parse", "--abbrev-ref", "HEAD"]),
            "ts_utc": datetime.now(UTC).isoformat(),
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
                "artifacts": record.artifacts,
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
    success_rate_ratio = summary["success_rate"]
    lines = [
        "# Test Summary",
        "",
        f"- Repo: `{meta['repo']}`",
        f"- Commit: {meta['commit']}",
        f"- Branch: {meta['branch']}",
        f"- Timestamp (UTC): {meta['ts_utc']}",
        f"- Duration: {meta['duration_sec']} s",
        f"- Success rate: {success_rate_ratio * 100:.2f}%",
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

    artifact_entries = [
        (test["nodeid"], test["artifacts"])
        for test in data["tests"]
        if test.get("artifacts")
    ]
    if artifact_entries:
        lines.append("")
        lines.append("## Test artifacts")
        for nodeid, artifacts in artifact_entries:
            for name, value in sorted(artifacts.items()):
                lines.append(f"- `{nodeid}`: {name} = `{value}`")
    else:
        lines.append("")
        lines.append("## Test artifacts")
        lines.append("- none")

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


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
    parser = argparse.ArgumentParser(description="Run pytest and emit reports")
    parser.add_argument(
        "--pytest-args",
        nargs=argparse.REMAINDER,
        help="Additional arguments forwarded to pytest",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging for the CLI and pytest log capture.",
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

    warnings.warn(DEPRECATION_MESSAGE, DeprecationWarning, stacklevel=2)
    logger.warning(DEPRECATION_MESSAGE)

    level = "DEBUG" if args.verbose else "INFO"
    base_logger_cfg = LoggerConfig(level=level, logger_name="run_tests")

    collector = ReportCollector()

    with setup_cli_logging("run_tests", base_logger_cfg) as logging_ctx:
        configure_logger(logging_ctx.log_cfg)

        extra_args = _normalise_extra_args(args.pytest_args)

        pytest_args = [
            "-q",
            f"--log-file={logging_ctx.log_path}",
            f"--log-file-level={logging_ctx.log_cfg.level}",
        ]

        if not _has_explicit_targets(extra_args):
            pytest_args.extend(_default_test_targets())

        pytest_args.extend(extra_args)

        pytest_exit_code = pytest.main(pytest_args, plugins=[collector])
        exit_code = int(pytest_exit_code)

        data = _build_json(collector, report_path=args.json)
        _write_markdown(args.markdown, data)

        success_rate_ratio = data["summary"]["success_rate"]
        if success_rate_ratio < SUCCESS_RATE_THRESHOLD:
            logger.error(
                "Success rate %.2f%% is below the required threshold of %.2f%%",
                success_rate_ratio * 100,
                SUCCESS_RATE_THRESHOLD * 100,
            )
            if exit_code == 0:
                exit_code = 1

    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
