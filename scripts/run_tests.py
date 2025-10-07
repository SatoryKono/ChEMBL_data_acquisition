"""Run the test suite and produce structured JSON and Markdown reports."""
from __future__ import annotations

import argparse
import json
import logging
import platform
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import pytest

from library.cli import configure_logger, create_logger_config
from library.cli.logging import setup_cli_logging

ROOT_DIR = Path(__file__).resolve().parents[1]
REPORTS_DIR = ROOT_DIR / "reports"
RAW_REPORT_FILE = REPORTS_DIR / "pytest_raw_report.json"
REPORT_FILE = REPORTS_DIR / "test_report.json"
SUMMARY_FILE = REPORTS_DIR / "test_summary.md"
COVERAGE_DIR = REPORTS_DIR / "coverage"
COVERAGE_XML = COVERAGE_DIR / "coverage.xml"
COVERAGE_HTML = COVERAGE_DIR / "html"
REPO_SLUG = "SatoryKono/ChEMBL_data_acquisition"
TEST_DIRECTORIES = (
    ROOT_DIR / "tests" / "unit",
    ROOT_DIR / "tests" / "integration",
    ROOT_DIR / "tests" / "postprocessing",
    ROOT_DIR / "tests" / "e2e",
)
_BASE_PYTEST_COMMAND: list[str] = [
    sys.executable,
    "-m",
    "pytest",
    "--json-report",
    "--json-report-file",
    str(RAW_REPORT_FILE),
    "--durations=0",
    "--cov=library",
    "--cov=scripts",
    "--cov-report=term",
    f"--cov-report=xml:{COVERAGE_XML}",
    f"--cov-report=html:{COVERAGE_HTML}",
    "-vv",
]
_DEFAULT_TEST_TARGETS: tuple[str, ...] = tuple(
    str(path) for path in TEST_DIRECTORIES if path.exists()
)


logger = logging.getLogger("run_tests")


def _relative_to_root(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT_DIR))
    except ValueError:
        return str(path)


def ensure_reports_directory() -> None:
    """Ensure the reports directory exists."""

    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    COVERAGE_DIR.mkdir(parents=True, exist_ok=True)


def run_pytest(command: Sequence[str]) -> int:
    """Execute ``command`` and stream output through the configured logger."""

    logger.debug("Executing pytest command: %s", " ".join(shlex.quote(part) for part in command))

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    assert process.stdout is not None  # for mypy; Popen(..., stdout=PIPE)
    with process.stdout:
        for raw_line in process.stdout:
            line = raw_line.rstrip()
            if line:
                logger.info(line)
            else:
                logger.info("")

    return_code = process.wait()
    if return_code != 0:
        logger.error("Pytest exited with non-zero status: %s", return_code)
    else:
        logger.debug("Pytest exited successfully")
    return return_code


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

    failure_rows = []
    for test in tests:
        status = str(test.get("status", "")).lower()
        if status not in {"failed", "error"}:
            continue
        nodeid = str(test.get("nodeid", "<unknown>"))
        message = _normalise_message(test.get("error"))
        failure_rows.append((nodeid, status, message))

    if not failure_rows:
        lines.append("- None")
    else:
        for nodeid, status, message in failure_rows:
            lines.append(f"- `{nodeid}` ({status})")
            display_message = message or "<no message>"
            lines.append("  ```")
            lines.extend(f"  {line}" for line in display_message.splitlines())
            lines.append("  ```")

    lines.append("")
    return "\n".join(lines)


def write_summary(report: dict[str, Any]) -> None:
    SUMMARY_FILE.write_text(build_summary_markdown(report), encoding="utf-8")


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the test suite and emit structured reports."
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Base logging level (default: INFO)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Shortcut for --log-level=DEBUG",
    )
    parser.add_argument(
        "--date",
        default=None,
        help="Date token forwarded to the log file suffix (format: YYYYMMDD)",
    )
    parser.add_argument(
        "pytest_args",
        nargs=argparse.REMAINDER,
        help="Arguments forwarded to pytest (use '--' before them)",
    )
    return parser.parse_args(argv)


def _extract_pytest_args(args: Sequence[str] | None) -> list[str]:
    if not args:
        return []
    if args and args[0] == "--":
        return list(args[1:])
    return list(args)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    level = str(args.log_level or "INFO").upper()
    if args.verbose:
        level = "DEBUG"

    log_cfg = create_logger_config(level)

    with setup_cli_logging("run_tests", log_cfg, args.date) as logging_ctx:
        configure_logger(logging_ctx.log_cfg)

        ensure_reports_directory()

        pytest_command = list(_BASE_PYTEST_COMMAND)
        pytest_command.extend(
            [
                "--log-file",
                str(logging_ctx.log_path),
                "--log-file-level",
                logging_ctx.log_cfg.level,
            ]
        )
        pytest_command.extend(_DEFAULT_TEST_TARGETS)
        pytest_command.extend(_extract_pytest_args(args.pytest_args))

        exit_code = run_pytest(pytest_command)

        raw_report = _load_raw_report()
        structured = build_structured_report(raw_report, exit_code)
        write_json_report(structured)
        write_summary(structured)

        log_path = _relative_to_root(logging_ctx.log_path)
        logger.info("Pytest finished with exit code %s", exit_code)
        logger.info("Log saved to %s", log_path)
        if RAW_REPORT_FILE.exists():
            logger.info(
                "Raw report available at %s",
                _relative_to_root(RAW_REPORT_FILE),
            )
        logger.info(
            "Structured report written to %s",
            _relative_to_root(REPORT_FILE),
        )
        logger.info(
            "Summary written to %s",
            _relative_to_root(SUMMARY_FILE),
        )

    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
