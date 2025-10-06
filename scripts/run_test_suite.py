"""Run the full pytest suite and export structured reports."""

from __future__ import annotations

import argparse
import json
import logging
import os
import platform
import subprocess
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import pytest


SUCCESS_RATE_THRESHOLD = 0.95
REPO_FALLBACK = "SatoryKono/ChEMBL_data_acquisition"


logger = logging.getLogger(__name__)


_STATUS_PRIORITY: dict[str, int] = {
    "passed": 0,
    "xpassed": 1,
    "skipped": 2,
    "xfailed": 3,
    "failed": 4,
    "error": 5,
    "pending": -1,
}


@dataclass
class TestResult:
    """Serializable representation of a single pytest outcome."""

    nodeid: str
    status: str = "pending"
    duration_ms: float = 0.0
    stdout: str = ""
    stderr: str = ""
    log: list[str] = field(default_factory=list)
    error: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "nodeid": self.nodeid,
            "status": self.status if self.status != "pending" else "error",
            "duration_ms": round(self.duration_ms, 3),
            "stdout": self.stdout,
            "stderr": self.stderr,
            "log": self.log,
            "error": self.error,
        }


class JsonReportPlugin:
    """Collect pytest results and persist optional per-suite logs."""

    def __init__(self, log_path: Path) -> None:
        self._log_path = log_path
        self._results: dict[str, TestResult] = {}
        self._started_at: float | None = None
        self._finished_at: float | None = None

    def pytest_sessionstart(self, session: pytest.Session) -> None:  # type: ignore[override]
        self._started_at = time.perf_counter()

    def pytest_sessionfinish(self, session: pytest.Session) -> None:  # type: ignore[override]
        self._finished_at = time.perf_counter()

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:  # type: ignore[override]
        if report.when not in {"setup", "call", "teardown"}:
            return

        entry = self._results.get(report.nodeid)
        if entry is None:
            entry = TestResult(nodeid=report.nodeid)
            self._results[report.nodeid] = entry

        entry.duration_ms += max(report.duration * 1000.0, 0.0)
        entry.stdout += getattr(report, "capstdout", "") or ""
        entry.stderr += getattr(report, "capstderr", "") or ""

        caplog = getattr(report, "caplog", None)
        if caplog:
            entry.log.append(caplog)

        for section_name, content in getattr(report, "sections", []):
            if "Captured log" in section_name:
                entry.log.append(content)

        status, message = _status_from_report(report)
        if status is None:
            return

        if _STATUS_PRIORITY.get(status, 0) >= _STATUS_PRIORITY.get(entry.status, -1):
            entry.status = status
            entry.error = message

    @property
    def results(self) -> list[TestResult]:
        ordered = sorted(self._results.values(), key=lambda item: item.nodeid)
        return ordered

    @property
    def duration_seconds(self) -> float:
        if self._started_at is None:
            return 0.0
        end = self._finished_at if self._finished_at is not None else time.perf_counter()
        return max(end - self._started_at, 0.0)


def _status_from_report(report: pytest.TestReport) -> tuple[str | None, str | None]:
    message: str | None = None
    status: str | None = None

    if report.when == "setup" and report.outcome == "failed":
        status = "error"
        message = _format_longrepr(report.longrepr)
    elif report.when == "teardown" and report.outcome == "failed":
        status = "error"
        message = _format_longrepr(report.longrepr)
    elif report.outcome == "skipped":
        status = "xfailed" if getattr(report, "wasxfail", None) else "skipped"
        message = _format_longrepr(report.longrepr)
    elif report.when == "call":
        if report.outcome == "failed":
            if getattr(report, "wasxfail", None):
                status = "xfailed"
                message = str(report.wasxfail)
            else:
                status = "failed"
                message = _format_longrepr(report.longrepr)
        elif report.outcome == "passed":
            if getattr(report, "wasxfail", None):
                status = "xpassed"
                message = str(report.wasxfail)
            else:
                status = "passed"
    return status, message


def _format_longrepr(longrepr: Any) -> str:
    if longrepr is None:
        return ""
    if isinstance(longrepr, tuple) and len(longrepr) == 3:
        return str(longrepr[2])
    if hasattr(longrepr, "longreprtext"):
        return str(longrepr.longreprtext)  # pragma: no cover - compatibility guard
    return str(longrepr)


def summarize_results(results: Sequence[TestResult]) -> dict[str, Any]:
    total = len(results)
    counts = {
        "passed": 0,
        "failed": 0,
        "skipped": 0,
        "xfailed": 0,
        "xpassed": 0,
        "error": 0,
    }
    for item in results:
        status = item.status if item.status != "pending" else "error"
        if status in counts:
            counts[status] += 1
    success_rate = 1.0 if total == 0 else counts["passed"] / total

    return {
        "total": total,
        **counts,
        "success_rate": success_rate,
    }


def _git_output(args: Sequence[str]) -> str:
    try:
        completed = subprocess.run(
            ["git", *args],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""
    return completed.stdout.strip()


def _normalise_repo_name(remote: str) -> str | None:
    remote = remote.strip()
    if not remote:
        return None
    if remote.endswith(".git"):
        remote = remote[:-4]
    if remote.startswith("git@"):
        parts = remote.split(":", 1)
        if len(parts) == 2:
            return parts[1]
    if remote.startswith("https://") or remote.startswith("http://"):
        parts = remote.split("//", 1)[-1].split("/", 1)
        if len(parts) == 2:
            return parts[1]
    if "/" in remote:
        return remote
    return None


def _gather_meta(duration_sec: float) -> dict[str, Any]:
    remote = _git_output(["config", "--get", "remote.origin.url"])
    repo_name = _normalise_repo_name(remote) or REPO_FALLBACK
    commit = _git_output(["rev-parse", "HEAD"]) or "unknown"
    branch = _git_output(["rev-parse", "--abbrev-ref", "HEAD"]) or "unknown"
    timestamp = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")

    return {
        "repo": repo_name,
        "commit": commit,
        "branch": branch,
        "ts_utc": timestamp,
        "duration_sec": round(duration_sec, 3),
        "python": platform.python_version(),
        "pytest": pytest.__version__,
    }


def _build_report(results: Sequence[TestResult], *, duration_sec: float) -> dict[str, Any]:
    summary = summarize_results(results)
    meta = _gather_meta(duration_sec)
    return {
        "meta": meta,
        "summary": summary,
        "tests": [result.to_dict() for result in results],
    }


def _write_json(report_path: Path, payload: dict[str, Any]) -> None:
    report_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _write_summary(summary_path: Path, payload: dict[str, Any]) -> None:
    meta = payload["meta"]
    summary = payload["summary"]

    success_rate = summary["success_rate"] * 100
    duration = meta.get("duration_sec", 0.0)

    lines = [
        "# Test Summary",
        "",
        f"- Repo: `{meta['repo']}`",
        f"- Commit: {meta['commit']}",
        f"- Branch: {meta['branch']}",
        f"- Timestamp (UTC): {meta['ts_utc']}",
        f"- Duration: {duration:.3f} s",
        f"- Success rate: {success_rate:.2f}%",
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
        "## Failed / Error details",
    ]

    failing = [
        (test["nodeid"], test.get("error") or test.get("stderr") or "")
        for test in payload["tests"]
        if test["status"] in {"failed", "error"}
    ]
    if failing:
        for nodeid, message in failing:
            short = message.splitlines()[0] if message else "(no message)"
            lines.append(f"- `{nodeid}`: `{short}`")
    else:
        lines.append("- none")

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
    pytest_exit_code = pytest.main(pytest_args, plugins=[plugin])
    exit_code = int(pytest_exit_code)

    results = plugin.results
    payload = _build_report(results, duration_sec=plugin.duration_seconds)
    summary = payload["summary"]

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
    _write_json(json_path, payload)
    _write_summary(summary_path, payload)

    return int(exit_code)


if __name__ == "__main__":
    raise SystemExit(main())
