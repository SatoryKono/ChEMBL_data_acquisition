"""Run the test suite and produce machine- and human-readable reports."""
from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

ROOT_DIR = Path(__file__).resolve().parents[1]
REPORTS_DIR = ROOT_DIR / "reports"
REPORT_FILE = REPORTS_DIR / "test_report.json"
LOG_FILE = REPORTS_DIR / "test_run.log"
SUMMARY_FILE = REPORTS_DIR / "test_summary.md"

PYTEST_COMMAND: tuple[str, ...] = (
    sys.executable,
    "-m",
    "pytest",
    "--json-report",
    "--json-report-file",
    str(REPORT_FILE),
    "--durations=0",
    "-vv",
)


def ensure_reports_directory() -> None:
    """Ensure that the reports directory exists."""

    REPORTS_DIR.mkdir(parents=True, exist_ok=True)


def run_pytest() -> int:
    """Execute pytest and capture its combined stdout/stderr into a log file."""

    with LOG_FILE.open("w", encoding="utf-8") as log_file:
        result = subprocess.run(
            PYTEST_COMMAND,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return result.returncode


def load_report() -> dict[str, Any]:
    """Load the JSON report produced by pytest-json-report, if available."""

    if not REPORT_FILE.exists():
        return {}
    try:
        return json.loads(REPORT_FILE.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _normalise_message(raw: Any) -> str:
    """Convert values from the report into a compact string."""

    if raw is None:
        return ""
    if isinstance(raw, str):
        return raw.strip()
    if isinstance(raw, Iterable) and not isinstance(raw, (bytes, bytearray)):
        joined = "\n".join(str(part) for part in raw)
        return joined.strip()
    return str(raw).strip()


def _extract_section_message(section: dict[str, Any]) -> str:
    """Extract a meaningful message from a report section."""

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


def collect_notable_failures(report: dict[str, Any], limit: int = 5) -> list[dict[str, str]]:
    """Return a concise list of failing test cases for the summary."""

    tests = report.get("tests", [])
    if not isinstance(tests, list):
        return []

    failures: list[dict[str, str]] = []
    for test in tests:
        if not isinstance(test, dict):
            continue
        outcome = test.get("outcome")
        if outcome in {"passed", "skipped", "xfailed"}:
            continue
        nodeid = str(test.get("nodeid", "<unknown>"))
        message = ""
        for phase in ("setup", "call", "teardown"):
            section = test.get(phase)
            if isinstance(section, dict):
                message = _extract_section_message(section)
            if message:
                break
        if not message:
            message = _normalise_message(test.get("longrepr"))
        if not message:
            message = "See test_run.log for details."
        truncated = message if len(message) <= 500 else f"{message[:497]}..."
        failures.append({
            "nodeid": nodeid,
            "outcome": str(outcome),
            "message": truncated,
        })
        if len(failures) >= limit:
            break
    return failures


def build_summary(report: dict[str, Any], exit_code: int) -> str:
    """Create a Markdown summary for the test run."""

    summary = report.get("summary", {}) if isinstance(report, dict) else {}
    total = int(summary.get("total", 0)) if summary else 0
    passed = int(summary.get("passed", 0)) if summary else 0
    failed = int(summary.get("failed", 0)) if summary else 0
    errors = int(summary.get("error", 0)) if summary else 0
    skipped = int(summary.get("skipped", 0)) if summary else 0
    xfailed = int(summary.get("xfailed", 0)) if summary else 0
    xpassed = int(summary.get("xpassed", 0)) if summary else 0

    success_rate = (passed / total * 100) if total else 0.0

    timestamp = datetime.now(timezone.utc).isoformat()

    lines = [
        "# Test Run Summary",
        "",
        f"- Timestamp (UTC): {timestamp}",
        f"- Exit code: {exit_code}",
        f"- JSON report: `{REPORT_FILE.relative_to(ROOT_DIR)}`",
        f"- Log file: `{LOG_FILE.relative_to(ROOT_DIR)}`",
    ]

    if report:
        duration = report.get("duration")
        if duration is not None:
            lines.append(f"- Duration: {duration:.2f}s")
    lines.extend(
        [
            "",
            "## Overview",
            "",
            "| Metric | Value |",
            "| --- | --- |",
            f"| Total tests | {total} |",
            f"| Passed | {passed} |",
            f"| Failed | {failed} |",
            f"| Errors | {errors} |",
            f"| Skipped | {skipped} |",
            f"| XFailed | {xfailed} |",
            f"| XPassed | {xpassed} |",
            f"| Success rate | {success_rate:.1f}% |",
        ]
    )

    failures = collect_notable_failures(report)
    lines.append("")
    if failures:
        lines.append("## Notable failures")
        lines.append("")
        for item in failures:
            lines.append(f"- `{item['nodeid']}` ({item['outcome']})")
            sanitized_message = item['message'].replace('\n', ' ')
            lines.append(f"  - {sanitized_message}")
        if len(failures) == 5:
            lines.append("")
            lines.append("_Only the first five failures are shown. Consult the log for full output._")
    else:
        if exit_code == 0:
            lines.append("All tests passed.")
        else:
            lines.append("No failure details captured. Check the log file for diagnostics.")

    lines.append("")
    return "\n".join(lines)


def write_summary(report: dict[str, Any], exit_code: int) -> None:
    """Persist the Markdown summary to disk."""

    SUMMARY_FILE.write_text(build_summary(report, exit_code), encoding="utf-8")


def main() -> int:
    """Run pytest, create artefacts, and forward the exit status."""

    ensure_reports_directory()
    exit_code = run_pytest()
    report = load_report()
    write_summary(report, exit_code)
    print(f"Pytest finished with exit code {exit_code}.")
    print(f"Log saved to {LOG_FILE.relative_to(ROOT_DIR)}.")
    if REPORT_FILE.exists():
        print(f"JSON report available at {REPORT_FILE.relative_to(ROOT_DIR)}.")
    else:
        print("JSON report was not created. See log for details.")
    print(f"Summary written to {SUMMARY_FILE.relative_to(ROOT_DIR)}.")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
