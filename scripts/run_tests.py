"""Run the test suite and produce structured JSON and Markdown reports."""

from __future__ import annotations

import argparse
import json
import logging
import os
import platform
import shlex
import shutil
import subprocess
import sys
from collections.abc import Mapping, MutableMapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
import threading
from typing import Any, cast
from xml.etree import ElementTree
from uuid import NAMESPACE_URL, uuid5

import pytest

from library.cli import configure_logger, create_logger_config
from library.cli.logging import setup_cli_logging
from library.cli_utils import resolve_invocation
from library.reporting.test_summary import (
    DEFAULT_REPO_SLUG,
    build_summary_markdown,
)
from library.reporting.test_summary import (
    normalise_message as _normalise_message,
)

ROOT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_REPORTS_DIR = ROOT_DIR / "reports"
RAW_REPORT_FILE = DEFAULT_REPORTS_DIR / "pytest_raw_report.json"
DEFAULT_REPORT_FILE = DEFAULT_REPORTS_DIR / "test_report.json"
DEFAULT_SUMMARY_FILE = DEFAULT_REPORTS_DIR / "test_summary.md"
COVERAGE_DIR = DEFAULT_REPORTS_DIR / "coverage"
COVERAGE_XML = COVERAGE_DIR / "coverage.xml"
COVERAGE_HTML = COVERAGE_DIR / "html"
# Backwards-compatible aliases for tests and external callers.
REPORTS_DIR = DEFAULT_REPORTS_DIR
REPORT_FILE = DEFAULT_REPORT_FILE
SUMMARY_FILE = DEFAULT_SUMMARY_FILE
try:
    DEFAULT_JSON_ARG = str(DEFAULT_REPORT_FILE.relative_to(ROOT_DIR))
except ValueError:  # pragma: no cover - defensive fallback
    DEFAULT_JSON_ARG = str(DEFAULT_REPORT_FILE)
try:
    DEFAULT_MARKDOWN_ARG = str(DEFAULT_SUMMARY_FILE.relative_to(ROOT_DIR))
except ValueError:  # pragma: no cover - defensive fallback
    DEFAULT_MARKDOWN_ARG = str(DEFAULT_SUMMARY_FILE)
REPO_SLUG = DEFAULT_REPO_SLUG
QUALITY_THRESHOLD_PERCENT = 95.0
LINE_COVERAGE_THRESHOLD_PERCENT = 80.0
QUALITY_FAILURE_EXIT_CODE = 1
VALIDATION_FAILURE_EXIT_CODE = 11
TEST_DIRECTORIES = (
    ROOT_DIR / "tests" / "unit",
    ROOT_DIR / "tests" / "integration",
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

PYTEST_TIMEOUT_ENV_VAR = "CHEMBL_DA_PYTEST_TIMEOUT"
_PROCESS_TERMINATION_TIMEOUT = 5.0


logger = logging.getLogger("run_tests")


class CoverageReportError(RuntimeError):
    """Raised when the coverage report cannot be parsed."""


def _relative_to_root(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT_DIR))
    except ValueError:
        return str(path)


def _clear_directory(directory: Path) -> None:
    """Remove all files and subdirectories inside ``directory``."""

    if not directory.exists():
        return

    for entry in directory.iterdir():
        try:
            if entry.is_dir():
                shutil.rmtree(entry)
            else:
                entry.unlink()
        except FileNotFoundError:  # pragma: no cover - race condition safety
            continue


def ensure_output_directories(report_file: Path, summary_file: Path) -> None:
    """Ensure directories for structured outputs exist."""

    REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    if COVERAGE_DIR.exists():
        _clear_directory(COVERAGE_DIR)
    COVERAGE_DIR.mkdir(parents=True, exist_ok=True)

    if COVERAGE_HTML.exists():
        _clear_directory(COVERAGE_HTML)
    COVERAGE_HTML.mkdir(parents=True, exist_ok=True)

    for path in (report_file, summary_file, RAW_REPORT_FILE):
        parent = path.parent
        parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            try:
                path.unlink()
            except IsADirectoryError:  # pragma: no cover - defensive fallback
                # A defensive guard in case an unexpected directory is present
                # at the target path. Removing directories is out of scope for
                # this helper, but we still make sure the call does not crash.
                pass


def _stream_process_output(process: subprocess.Popen[str]) -> None:
    """Stream stdout from ``process`` to the configured logger."""

    stdout = process.stdout
    if stdout is None:  # pragma: no cover - defensive fallback
        return
    with stdout:
        for raw_line in stdout:
            line = raw_line.rstrip()
            if line:
                logger.info(line)
            else:
                logger.info("")


def run_pytest(command: Sequence[str], *, timeout: float | None = None) -> int:
    """Execute ``command`` and stream output through the configured logger."""

    logger.debug(
        "Executing pytest command: %s", " ".join(shlex.quote(part) for part in command)
    )

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    stream_thread = threading.Thread(
        target=_stream_process_output, args=(process,), daemon=True
    )
    stream_thread.start()

    timed_out = False
    try:
        return_code = process.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        logger.error(
            "Pytest command exceeded timeout of %s seconds; terminating process",
            timeout,
        )
        process.kill()
        try:
            return_code = process.wait(timeout=_PROCESS_TERMINATION_TIMEOUT)
        except subprocess.TimeoutExpired:  # pragma: no cover - defensive fallback
            logger.error(
                "Pytest process did not exit after kill; waiting without timeout"
            )
            return_code = process.wait()
    finally:
        stream_thread.join(timeout=_PROCESS_TERMINATION_TIMEOUT)
        if stream_thread.is_alive():  # pragma: no cover - defensive fallback
            logger.debug("Pytest output stream thread is still running after join")

    if timed_out:
        logger.error("Pytest run aborted due to timeout")
    if return_code != 0:
        logger.error("Pytest exited with non-zero status: %s", return_code)
    else:
        logger.debug("Pytest exited successfully")
    return return_code


def _parse_line_coverage(path: Path) -> float:
    """Return the total line coverage percentage from ``coverage.xml``."""

    if not path.exists():
        raise CoverageReportError(
            f"Coverage report {_relative_to_root(path)} was not generated"
        )

    try:
        tree = ElementTree.parse(path)
    except ElementTree.ParseError as exc:  # pragma: no cover - depends on coverage output
        raise CoverageReportError(
            f"Coverage report {_relative_to_root(path)} is not valid XML: {exc}"
        ) from exc

    root = tree.getroot()
    line_rate_raw = root.attrib.get("line-rate")
    if line_rate_raw is None:
        raise CoverageReportError(
            "Coverage XML is missing the 'line-rate' attribute on the root element"
        )

    try:
        line_rate = float(line_rate_raw)
    except ValueError as exc:
        raise CoverageReportError(
            "Coverage XML attribute 'line-rate' must be numeric"
        ) from exc

    if line_rate > 1.0:
        return max(0.0, min(100.0, line_rate))
    return max(0.0, min(100.0, line_rate * 100.0))


def _load_raw_report() -> dict[str, Any]:
    if not RAW_REPORT_FILE.exists():
        return {}
    try:
        return json.loads(RAW_REPORT_FILE.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
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


def _calculate_success_rate(summary: Mapping[str, int | float]) -> float:
    total = int(summary.get("total", 0) or 0)
    skipped = int(summary.get("skipped", 0) or 0)
    passed = int(summary.get("passed", 0) or 0)
    xfailed = int(summary.get("xfailed", 0) or 0)

    executed_total = max(0, total - skipped)
    if executed_total == 0:
        return 0.0

    successes = passed + xfailed
    ratio = successes / executed_total
    return max(0.0, min(1.0, ratio))


def build_structured_report(raw: dict[str, Any], exit_code: int) -> dict[str, Any]:
    tests_raw = raw.get("tests", [])
    tests: list[dict[str, Any]] = []
    summary: dict[str, int | float] = {
        "total": 0,
        "passed": 0,
        "failed": 0,
        "skipped": 0,
        "xfailed": 0,
        "xpassed": 0,
        "error": 0,
        "success_rate": 0.0,
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

    success_rate = _calculate_success_rate(summary)

    if summary["total"] == 0 and exit_code != 0:
        message = (
            "pytest exited with code" f" {exit_code} before completing the test run"
        )
        tests.append(
            {
                "nodeid": "<pytest-startup>",
                "status": "error",
                "duration_ms": 0.0,
                "stdout": "",
                "stderr": "",
                "log": [],
                "error": message,
            }
        )
        summary["total"] = 1
        summary["passed"] = 0
        summary["failed"] = 0
        summary["skipped"] = 0
        summary["xfailed"] = 0
        summary["xpassed"] = 0
        summary["error"] = 1
        summary["success_rate"] = 0.0
    else:
        summary["success_rate"] = round(success_rate, 4)

    meta = {
        "repo": REPO_SLUG,
        "commit": _git_output("rev-parse", "HEAD"),
        "branch": _git_output("rev-parse", "--abbrev-ref", "HEAD"),
        "ts_utc": datetime.now(UTC).isoformat(),
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


def _ensure_failure_record(
    report: dict[str, Any], exit_code: int, reason: str
) -> dict[str, Any]:
    summary = cast(MutableMapping[str, object], report.setdefault("summary", {}))
    raw_total = summary.get("total", 0)
    total_value = int(raw_total) if isinstance(raw_total, (int, float)) else 0
    summary["total"] = max(total_value, 1)
    summary["passed"] = 0
    summary["failed"] = 0
    summary["skipped"] = 0
    summary["xfailed"] = 0
    summary["xpassed"] = 0
    raw_error = summary.get("error", 0)
    error_value = int(raw_error) if isinstance(raw_error, (int, float)) else 0
    summary["error"] = max(error_value, 1)
    summary["success_rate"] = 0.0

    entry = {
        "nodeid": "<report-generation>",
        "status": "error",
        "duration_ms": 0.0,
        "stdout": "",
        "stderr": "",
        "log": [],
        "error": reason or f"failed to build structured report (exit code {exit_code})",
    }

    tests = report.setdefault("tests", [])
    if entry not in tests:
        tests.append(entry)
    meta = report.setdefault("meta", {})
    meta.setdefault("exit_code", exit_code)
    meta.setdefault("repo", REPO_SLUG)
    meta.setdefault("branch", _git_output("rev-parse", "--abbrev-ref", "HEAD"))
    meta.setdefault("commit", _git_output("rev-parse", "HEAD"))
    meta.setdefault("ts_utc", datetime.now(UTC).isoformat())
    meta.setdefault("duration_sec", 0.0)
    meta.setdefault("python", platform.python_version())
    meta.setdefault("pytest", pytest.__version__)
    return report


def write_json_report(report: dict[str, Any], destination: Path) -> None:
    destination.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def validate_structured_report(report: dict[str, Any]) -> None:
    if not isinstance(report, dict):
        raise ValueError("Report must be a dictionary")

    summary = report.get("summary")
    if not isinstance(summary, dict):
        raise ValueError("Report summary is missing or malformed")

    required_summary_keys = {
        "total",
        "passed",
        "failed",
        "skipped",
        "xfailed",
        "xpassed",
        "error",
        "success_rate",
    }
    missing_keys = required_summary_keys - set(summary)
    if missing_keys:
        raise ValueError(f"Report summary missing keys: {sorted(missing_keys)}")

    for key in required_summary_keys - {"success_rate"}:
        value = summary.get(key)
        if not isinstance(value, int) or value < 0:
            raise ValueError(f"Summary field '{key}' must be a non-negative integer")

    success_rate = summary.get("success_rate")
    if not isinstance(success_rate, int | float):
        raise ValueError("Summary field 'success_rate' must be numeric")
    if not 0.0 <= float(success_rate) <= 1.0:
        raise ValueError("Summary field 'success_rate' must be between 0 and 1")

    tests = report.get("tests")
    if not isinstance(tests, list):
        raise ValueError("Report tests section must be a list")
    for index, entry in enumerate(tests):
        if not isinstance(entry, dict):
            raise ValueError(f"Test entry at index {index} must be a dictionary")
        if "nodeid" not in entry or not isinstance(entry["nodeid"], str):
            raise ValueError(f"Test entry at index {index} missing string 'nodeid'")
        if "status" not in entry or not isinstance(entry["status"], str):
            raise ValueError(f"Test entry at index {index} missing string 'status'")
        duration = entry.get("duration_ms")
        if not isinstance(duration, int | float) or duration < 0:
            raise ValueError(
                f"Test entry at index {index} has invalid 'duration_ms' (must be >= 0)"
            )
        for text_field in ("stdout", "stderr"):
            if text_field in entry and not isinstance(entry[text_field], str):
                raise ValueError(
                    f"Test entry at index {index} has non-string '{text_field}'"
                )
        log_entries = entry.get("log", [])
        if not isinstance(log_entries, list) or not all(
            isinstance(log_entry, str) for log_entry in log_entries
        ):
            raise ValueError(
                f"Test entry at index {index} must contain a list of string logs"
            )
        error_field = entry.get("error")
        if error_field is not None and not isinstance(error_field, str):
            raise ValueError(
                f"Test entry at index {index} has non-string 'error' field"
            )


def validate_report_file(path: Path) -> None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ValueError(f"Report file {path} was not created") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"Report file {path} contains invalid JSON") from exc
    validate_structured_report(payload)
def write_summary(report: dict[str, Any], destination: Path) -> None:
    destination.write_text(build_summary_markdown(report), encoding="utf-8")


def _log_run_artifacts(
    logging_ctx: Any, exit_code: int, report_path: Path, summary_path: Path
) -> None:
    """Emit log lines describing where artefacts were written."""

    logger.info("Pytest finished with exit code %s", exit_code)

    log_path = getattr(logging_ctx, "log_path", None)
    if isinstance(log_path, Path):
        logger.info("Log saved to %s", _relative_to_root(log_path))
    elif isinstance(log_path, str):
        logger.info("Log saved to %s", log_path)

    if RAW_REPORT_FILE.exists():
        logger.info("Raw report available at %s", _relative_to_root(RAW_REPORT_FILE))
    if report_path.exists():
        logger.info(
            "Structured report written to %s",
            _relative_to_root(report_path),
        )
    if summary_path.exists():
        logger.info(
            "Summary written to %s",
            _relative_to_root(summary_path),
        )


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
        "--json",
        dest="json_path",
        default=DEFAULT_JSON_ARG,
        help="Path to the structured JSON report (relative paths resolve from repo root)",
    )
    parser.add_argument(
        "--markdown",
        dest="markdown_path",
        default=DEFAULT_MARKDOWN_ARG,
        help="Path to the Markdown summary report (relative paths resolve from repo root)",
    )
    parser.add_argument(
        "--run-id",
        dest="run_id",
        default=os.environ.get("CHEMBL_DA_RUN_ID"),
        help="Override the run identifier used for logging",
    )
    parser.add_argument(
        "--pytest-timeout",
        dest="pytest_timeout",
        type=float,
        default=None,
        help=(
            "Maximum number of seconds to allow pytest to run before forcing termination. "
            "Can also be provided via the CHEMBL_DA_PYTEST_TIMEOUT environment variable."
        ),
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


def _resolve_output_path(raw: str | None, default: Path) -> Path:
    if not raw:
        return default
    path = Path(raw)
    if not path.is_absolute():
        path = ROOT_DIR / path
    return path


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    invocation = resolve_invocation("run_tests", argv)
    level = str(args.log_level or "INFO").upper()
    if args.verbose:
        level = "DEBUG"

    report_path = _resolve_output_path(args.json_path, DEFAULT_REPORT_FILE)
    summary_path = _resolve_output_path(args.markdown_path, DEFAULT_SUMMARY_FILE)

    run_id_value = args.run_id.strip() if isinstance(args.run_id, str) else args.run_id
    if isinstance(run_id_value, str) and not run_id_value:
        run_id_value = None
    if not run_id_value:
        descriptor = "\n".join(
            [
                *invocation,
                f"json={report_path.resolve()}",
                f"markdown={summary_path.resolve()}",
            ]
        )
        run_id_value = uuid5(NAMESPACE_URL, descriptor).hex

    log_cfg = create_logger_config(level, run_id=run_id_value)

    with setup_cli_logging("run_tests", log_cfg, args.date) as logging_ctx:
        configure_logger(logging_ctx.log_cfg)

        ensure_output_directories(report_path, summary_path)

        timeout_value = args.pytest_timeout
        if timeout_value is None:
            env_timeout = os.environ.get(PYTEST_TIMEOUT_ENV_VAR)
            if env_timeout:
                try:
                    timeout_value = float(env_timeout)
                except ValueError:
                    logger.warning(
                        "Ignoring invalid %s value: %r",
                        PYTEST_TIMEOUT_ENV_VAR,
                        env_timeout,
                    )
        if timeout_value is not None and timeout_value <= 0:
            logger.warning(
                "Ignoring non-positive pytest timeout value: %s",
                timeout_value,
            )
            timeout_value = None
        if timeout_value is not None:
            logger.info("Pytest timeout configured to %.2f seconds", timeout_value)

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

        exit_code = run_pytest(pytest_command, timeout=timeout_value)

        structured: dict[str, Any] | None = None
        validation_exit_code: int | None = None
        generation_error: str | None = None

        try:
            raw_report = _load_raw_report()
            structured = build_structured_report(raw_report, exit_code)
            validate_structured_report(structured)
        except ValueError as exc:  # pragma: no cover - defensive guard
            logger.error("Structured report validation failed: %s", exc)
            validation_exit_code = VALIDATION_FAILURE_EXIT_CODE
            generation_error = _normalise_message(exc)
        except Exception as exc:  # pragma: no cover - defensive guard
            logger.error("Failed to build structured report: %s", exc)
            validation_exit_code = VALIDATION_FAILURE_EXIT_CODE
            generation_error = _normalise_message(exc)
        finally:
            if structured is None:
                structured = build_structured_report({}, exit_code)
            if generation_error:
                structured = _ensure_failure_record(
                    structured,
                    exit_code,
                    generation_error,
                )

            try:
                write_json_report(structured, report_path)
            except Exception as exc:  # pragma: no cover - defensive guard
                logger.error(
                    "Failed to write structured report to %s: %s",
                    _relative_to_root(report_path),
                    exc,
                )
            else:
                try:
                    validate_report_file(report_path)
                except ValueError as exc:  # pragma: no cover - defensive guard
                    logger.error("Written report failed validation: %s", exc)
                    validation_exit_code = VALIDATION_FAILURE_EXIT_CODE
        try:
            write_summary(structured, summary_path)
        except Exception as exc:  # pragma: no cover - defensive guard
            logger.error(
                "Failed to write summary to %s: %s",
                _relative_to_root(summary_path),
                exc,
            )
            _log_run_artifacts(logging_ctx, exit_code, report_path, summary_path)
            return VALIDATION_FAILURE_EXIT_CODE

        _log_run_artifacts(logging_ctx, exit_code, report_path, summary_path)

        final_exit_code = exit_code
        success_rate_raw = structured.get("summary", {}).get("success_rate", 0.0) or 0.0
        try:
            success_rate_value = float(success_rate_raw)
        except (TypeError, ValueError):  # pragma: no cover - guarded by validation
            logger.error(
                "Structured summary provided a non-numeric success rate %r; treating it as 0%%",
                success_rate_raw,
            )
            success_rate_value = 0.0

        success_rate_pct = (
            success_rate_value * 100.0
            if success_rate_value <= 1.0
            else success_rate_value
        )

        if success_rate_pct < QUALITY_THRESHOLD_PERCENT:
            logger.error(
                "Success rate %.2f%% is below the required %.2f%% threshold",
                success_rate_pct,
                QUALITY_THRESHOLD_PERCENT,
            )
            if final_exit_code == 0:
                final_exit_code = QUALITY_FAILURE_EXIT_CODE
        else:
            logger.info(
                "Success rate %.2f%% meets the required %.2f%% threshold",
                success_rate_pct,
                QUALITY_THRESHOLD_PERCENT,
            )

        if validation_exit_code is not None:
            final_exit_code = validation_exit_code

        coverage_exit_code: int | None = None
        try:
            coverage_pct = _parse_line_coverage(COVERAGE_XML)
        except CoverageReportError as exc:
            logger.error("Unable to determine coverage: %s", exc)
            coverage_exit_code = QUALITY_FAILURE_EXIT_CODE
        else:
            if coverage_pct < LINE_COVERAGE_THRESHOLD_PERCENT:
                logger.error(
                    "Line coverage %.2f%% is below the required %.2f%% threshold",
                    coverage_pct,
                    LINE_COVERAGE_THRESHOLD_PERCENT,
                )
                coverage_exit_code = QUALITY_FAILURE_EXIT_CODE
            else:
                logger.info(
                    "Line coverage %.2f%% meets the required %.2f%% threshold",
                    coverage_pct,
                    LINE_COVERAGE_THRESHOLD_PERCENT,
                )

        if coverage_exit_code is not None and final_exit_code == 0:
            final_exit_code = coverage_exit_code

    return final_exit_code


if __name__ == "__main__":
    raise SystemExit(main())
