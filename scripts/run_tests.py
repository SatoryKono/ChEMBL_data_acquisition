"""Run the test suite and produce structured JSON and Markdown reports."""
from __future__ import annotations

import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pytest

ROOT_DIR = Path(__file__).resolve().parents[1]
REPORTS_DIR = ROOT_DIR / "reports"
RAW_REPORT_FILE = REPORTS_DIR / "pytest_raw_report.json"
REPORT_FILE = REPORTS_DIR / "test_report.json"
LOG_FILE = REPORTS_DIR / "test_run.log"
SUMMARY_FILE = REPORTS_DIR / "test_summary.md"
REPO_SLUG = "SatoryKono/ChEMBL_data_acquisition"

PYTEST_COMMAND: tuple[str, ...] = (
    sys.executable,
    "-m",
    "pytest",
    "--json-report",
    "--json-report-file",
    str(RAW_REPORT_FILE),
    "--durations=0",
    "-vv",
)


def ensure_reports_directory() -> None:
    """Ensure the reports directory exists."""

    REPORTS_DIR.mkdir(parents=True, exist_ok=True)


def run_pytest() -> int:
    """Execute pytest and capture combined output in the log file."""

    with LOG_FILE.open("w", encoding="utf-8") as log_file:
        result = subprocess.run(
            PYTEST_COMMAND,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return result.returncode


def _load_raw_report() -> dict[str, Any]:
    if not RAW_REPORT_FILE.exists():
        return {}
    try:
        return json.loads(RAW_REPORT_FILE.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _normalise_message(raw: Any) -> str:
    if raw is None:
        return ""
    if isinstance(raw, str):
        return raw.strip()
    if isinstance(raw, Iterable) and not isinstance(raw, (bytes, bytearray)):
        joined = "\n".join(str(part) for part in raw)
        return joined.strip()
    return str(raw).strip()


def _extract_section_message(section: dict[str, Any]) -> str:
    for key in ("longrepr", "message"):
        if key in section:
            message = _normalise_message(section[key])
            if message:
                return message
    crash = section.get("crash")
    if isinstance(crash, dict):
        for key in ("message", "path"):
            value = crash.get(key)
            if value:
                crash_message = _normalise_message(value)
                if crash_message:
                    return crash_message
    return ""


def _git_output(*args: str) -> str:
    try:
        return subprocess.check_output(
            ("git", *args),
            cwd=ROOT_DIR,
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _collect_phase_output(section: dict[str, Any], field: str) -> list[str]:
    if not isinstance(section, dict):
        return []
    value = section.get(field)
    if value is None:
        return []
    if isinstance(value, list):
        return [_normalise_message(item) for item in value if _normalise_message(item)]
    text = _normalise_message(value)
    return [text] if text else []


def _build_test_entry(test: dict[str, Any]) -> dict[str, Any]:
    nodeid = str(test.get("nodeid", "<unknown>"))
    status = str(test.get("outcome", "unknown")).lower()
    duration = 0.0
    stdout_parts: list[str] = []
    stderr_parts: list[str] = []
    log_entries: list[str] = []

    for phase_name in ("setup", "call", "teardown"):
        section = test.get(phase_name)
        if not isinstance(section, dict):
            continue
        duration += float(section.get("duration", 0.0) or 0.0)
        stdout_parts.extend(_collect_phase_output(section, "stdout"))
        stderr_parts.extend(_collect_phase_output(section, "stderr"))
        log_entries.extend(_collect_phase_output(section, "log"))

    error_message = ""
    if status not in {"passed", "skipped", "xfailed", "xpassed"}:
        for phase_name in ("setup", "call", "teardown"):
            section = test.get(phase_name)
            if isinstance(section, dict):
                error_message = _extract_section_message(section)
            if error_message:
                break
        if not error_message:
            error_message = _normalise_message(test.get("longrepr"))

    entry = {
        "nodeid": nodeid,
        "status": status or "unknown",
        "duration_ms": round(duration * 1000.0, 3),
        "stdout": "\n".join(stdout_parts),
        "stderr": "\n".join(stderr_parts),
        "log": log_entries,
        "error": error_message or None,
    }
    return entry


def build_structured_report(raw: dict[str, Any], exit_code: int) -> dict[str, Any]:
    tests_raw = raw.get("tests", [])
    tests: list[dict[str, Any]] = []
    summary = {
        "total": 0,
        "passed": 0,
        "failed": 0,
        "skipped": 0,
        "xfailed": 0,
        "xpassed": 0,
        "error": 0,
    }

    if isinstance(tests_raw, list):
        iterable = tests_raw
    else:
        iterable = []

    for test in iterable:
        if not isinstance(test, dict):
            continue
        entry = _build_test_entry(test)
        tests.append(entry)
        summary["total"] += 1
        status = entry["status"]
        if status in summary:
            summary[status] += 1
        elif status.startswith("xfail"):
            summary["xfailed"] += 1
        elif status.startswith("xpass"):
            summary["xpassed"] += 1
        else:
            summary["error"] += 1

    success_rate = (
        (summary["passed"] / summary["total"] * 100.0)
        if summary["total"]
        else 0.0
    )
    summary["success_rate"] = round(success_rate, 2)

    meta = {
        "repo": REPO_SLUG,
        "commit": _git_output("rev-parse", "HEAD"),
        "branch": _git_output("rev-parse", "--abbrev-ref", "HEAD"),
        "ts_utc": datetime.now(timezone.utc).isoformat(),
        "duration_sec": float(raw.get("duration", 0.0) or 0.0),
        "python": platform.python_version(),
        "pytest": pytest.__version__,
        "exit_code": exit_code,
    }

    return {
        "meta": meta,
        "summary": summary,
        "tests": tests,
    }


def write_json_report(report: dict[str, Any]) -> None:
    REPORT_FILE.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def build_summary_markdown(report: dict[str, Any]) -> str:
    meta = report.get("meta", {})
    summary = report.get("summary", {})
    tests = report.get("tests", [])

    repo = meta.get("repo", REPO_SLUG)
    commit = meta.get("commit", "unknown")
    branch = meta.get("branch", "unknown")
    timestamp = meta.get("ts_utc", datetime.now(timezone.utc).isoformat())
    duration = float(meta.get("duration_sec", 0.0) or 0.0)
    success_rate = summary.get("success_rate", 0.0)

    lines = [
        "# Test Summary",
        "",
        f"- Repo: `{repo}`",
        f"- Commit: {commit}",
        f"- Branch: {branch}",
        f"- Timestamp (UTC): {timestamp}",
        f"- Duration: {duration:.2f} s",
        f"- Success rate: {success_rate:.2f}%",
        "",
        "| total | passed | failed | skipped | xfailed | xpassed | error |",
        "|------:|-------:|-------:|--------:|--------:|--------:|------:|",
        "| {total:5d} | {passed:5d} | {failed:5d} | {skipped:6d} | {xfailed:6d} | {xpassed:6d} | {error:5d} |".format(
            total=summary.get("total", 0),
            passed=summary.get("passed", 0),
            failed=summary.get("failed", 0),
            skipped=summary.get("skipped", 0),
            xfailed=summary.get("xfailed", 0),
            xpassed=summary.get("xpassed", 0),
            error=summary.get("error", 0),
        ),
        "",
        "## Failed / Error details",
    ]

    failure_rows = [
        (test["nodeid"], (test.get("error") or "See reports/test_run.log"))
        for test in tests
        if test.get("status") in {"failed", "error"}
    ]

    if not failure_rows:
        lines.append("- None")
    else:
        for nodeid, message in failure_rows:
            compact_message = message.replace("\n", " ").strip()
            lines.append(f"- `{nodeid}`: {compact_message}")

    lines.append("")
    return "\n".join(lines)


def write_summary(report: dict[str, Any]) -> None:
    SUMMARY_FILE.write_text(build_summary_markdown(report), encoding="utf-8")


def main() -> int:
    ensure_reports_directory()
    exit_code = run_pytest()
    raw_report = _load_raw_report()
    structured = build_structured_report(raw_report, exit_code)
    write_json_report(structured)
    write_summary(structured)

    print(f"Pytest finished with exit code {exit_code}.")
    print(f"Log saved to {LOG_FILE.relative_to(ROOT_DIR)}.")
    if RAW_REPORT_FILE.exists():
        print(f"Raw report available at {RAW_REPORT_FILE.relative_to(ROOT_DIR)}.")
    print(f"Structured report written to {REPORT_FILE.relative_to(ROOT_DIR)}.")
    print(f"Summary written to {SUMMARY_FILE.relative_to(ROOT_DIR)}.")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
