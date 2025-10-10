"""Unit tests for :mod:`scripts.run_tests`."""

from __future__ import annotations

import json
import logging
import sys
from collections.abc import Sequence
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from library.cli import LoggerConfig
from scripts import run_tests


@pytest.mark.unit
def test_run_tests__verbose_creates_debug_log(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reports_dir = tmp_path / "reports"
    coverage_dir = reports_dir / "coverage"
    raw_report_file = reports_dir / "pytest_raw_report.json"
    report_file = reports_dir / "test_report.json"
    summary_file = reports_dir / "test_summary.md"

    monkeypatch.setattr(run_tests, "ROOT_DIR", tmp_path)
    monkeypatch.setattr(run_tests, "REPORTS_DIR", reports_dir, raising=False)
    monkeypatch.setattr(run_tests, "RAW_REPORT_FILE", raw_report_file, raising=False)
    monkeypatch.setattr(run_tests, "REPORT_FILE", report_file, raising=False)
    monkeypatch.setattr(run_tests, "SUMMARY_FILE", summary_file, raising=False)
    monkeypatch.setattr(run_tests, "COVERAGE_DIR", coverage_dir, raising=False)
    monkeypatch.setattr(
        run_tests, "COVERAGE_XML", coverage_dir / "coverage.xml", raising=False
    )
    monkeypatch.setattr(
        run_tests, "COVERAGE_HTML", coverage_dir / "html", raising=False
    )

    base_command: list[str] = [
        sys.executable,
        "-m",
        "pytest",
        "--json-report",
        "--json-report-file",
        str(raw_report_file),
        "--durations=0",
        "--cov=library",
        "--cov=scripts",
        "--cov-report=term",
        f"--cov-report=xml:{coverage_dir / 'coverage.xml'}",
        f"--cov-report=html:{coverage_dir / 'html'}",
        "-vv",
    ]
    monkeypatch.setattr(run_tests, "_BASE_PYTEST_COMMAND", base_command, raising=False)
    monkeypatch.setattr(
        run_tests, "_DEFAULT_TEST_TARGETS", ("tests/unit",), raising=False
    )
    monkeypatch.setattr(run_tests, "_git_output", lambda *args, **kwargs: "stub")

    captured_log_path: Path | None = None

    @contextmanager
    def _fake_setup(script_name: str, log_cfg: LoggerConfig, date: str | None = None):
        nonlocal captured_log_path
        assert script_name == "run_tests"
        assert date == "20240131"
        captured_log_path = tmp_path / "logs" / f"run_tests_{date}.log"
        cloned_cfg = LoggerConfig(
            level=log_cfg.level,
            run_id=log_cfg.run_id,
            redact_secrets=log_cfg.redact_secrets,
            stream=log_cfg.stream,
            handlers=list(log_cfg.handlers),
            logger_name=log_cfg.logger_name,
        )
        yield SimpleNamespace(
            log_path=captured_log_path,
            log_cfg=cloned_cfg,
            console_stream=None,
        )

    captured_commands: list[list[str]] = []

    def _fake_run_pytest(command: Sequence[str]) -> int:
        captured_commands.append(list(command))
        return 0

    captured_configs: list[LoggerConfig] = []

    def _fake_configure(cfg: LoggerConfig) -> object:
        captured_configs.append(cfg)
        return object()

    monkeypatch.setattr(run_tests, "run_pytest", _fake_run_pytest)
    monkeypatch.setattr(run_tests, "configure_logger", _fake_configure)
    monkeypatch.setattr(run_tests, "setup_cli_logging", _fake_setup)
    monkeypatch.setattr(
        run_tests,
        "_load_raw_report",
        lambda: {"tests": [], "duration": 0.0},
    )

    exit_code = run_tests.main(
        [
            "--verbose",
            "--date",
            "20240131",
            "--",
            "-k",
            "unit",
        ]
    )

    assert exit_code == 0
    assert captured_commands, "run_pytest should be invoked"

    command = captured_commands[-1]
    assert "--log-file" in command
    log_path_token = command[command.index("--log-file") + 1]
    assert captured_log_path is not None
    expected_log_path = tmp_path / "logs" / "run_tests_20240131.log"
    assert captured_log_path == expected_log_path
    assert Path(log_path_token) == captured_log_path

    assert "--log-file-level" in command
    level_value = command[command.index("--log-file-level") + 1]
    assert level_value == "DEBUG"

    assert command[-2:] == ["-k", "unit"]

    assert report_file.exists()
    assert summary_file.exists()
    assert captured_configs and captured_configs[-1].level == "DEBUG"


@pytest.mark.unit
def test_build_structured_report__captures_failure_messages() -> None:
    raw_report = {
        "duration": 1.5,
        "tests": [
            {
                "nodeid": "tests/unit/test_example.py::test_failure",
                "outcome": "failed",
                "call": {
                    "duration": 0.25,
                    "longrepr": "AssertionError: boom\nline 1\nline 2",
                    "stdout": "",
                    "stderr": "",
                    "log": [],
                },
            },
            {
                "nodeid": "tests/unit/test_example.py::test_success",
                "outcome": "passed",
                "call": {
                    "duration": 0.75,
                    "stdout": "",
                    "stderr": "",
                    "log": [],
                },
            },
        ],
    }

    structured = run_tests.build_structured_report(raw_report, exit_code=1)

    assert structured["meta"]["exit_code"] == 1
    assert structured["summary"] == {
        "total": 2,
        "passed": 1,
        "failed": 1,
        "skipped": 0,
        "xfailed": 0,
        "xpassed": 0,
        "error": 0,
        "success_rate": pytest.approx(0.5),
    }

    failure_entry = next(
        item for item in structured["tests"] if item["nodeid"].endswith("test_failure")
    )
    assert failure_entry["status"] == "failed"
    assert failure_entry["error"] == "AssertionError: boom\nline 1\nline 2"


@pytest.mark.unit
def test_build_structured_report__injects_placeholder_on_startup_failure() -> None:
    raw_report: dict[str, Any] = {"duration": 0.0, "tests": []}

    structured = run_tests.build_structured_report(raw_report, exit_code=2)

    summary = structured["summary"]
    assert summary["total"] == 1
    assert summary["error"] == 1
    assert summary["success_rate"] == 0.0

    placeholder = structured["tests"][0]
    assert placeholder["status"] == "error"
    assert "code 2" in placeholder["error"]
@pytest.mark.unit
def test_calculate_success_rate__counts_passed_and_xfailed() -> None:
    summary = {
        "total": 6,
        "passed": 3,
        "failed": 1,
        "skipped": 1,
        "xfailed": 1,
        "xpassed": 0,
        "error": 0,
    }

    success_rate = run_tests._calculate_success_rate(summary)

    # Executed tests exclude the skipped item -> 5 executed, 4 successes
    assert success_rate == pytest.approx(4 / 5)


@pytest.mark.unit
def test_calculate_success_rate__all_skipped_return_unity() -> None:
    summary = {
        "total": 3,
        "passed": 0,
        "failed": 0,
        "skipped": 3,
        "xfailed": 0,
        "xpassed": 0,
        "error": 0,
    }

    success_rate = run_tests._calculate_success_rate(summary)

    assert success_rate == pytest.approx(1.0)


@pytest.mark.unit
def test_main__returns_exit_code_one_when_quality_gate_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    reports_dir = tmp_path / "reports"
    coverage_dir = reports_dir / "coverage"
    raw_report_file = reports_dir / "pytest_raw_report.json"
    report_file = reports_dir / "test_report.json"
    summary_file = reports_dir / "test_summary.md"

    monkeypatch.setattr(run_tests, "ROOT_DIR", tmp_path, raising=False)
    monkeypatch.setattr(run_tests, "REPORTS_DIR", reports_dir, raising=False)
    monkeypatch.setattr(run_tests, "RAW_REPORT_FILE", raw_report_file, raising=False)
    monkeypatch.setattr(run_tests, "REPORT_FILE", report_file, raising=False)
    monkeypatch.setattr(run_tests, "SUMMARY_FILE", summary_file, raising=False)
    monkeypatch.setattr(run_tests, "COVERAGE_DIR", coverage_dir, raising=False)
    monkeypatch.setattr(
        run_tests, "COVERAGE_XML", coverage_dir / "coverage.xml", raising=False
    )
    monkeypatch.setattr(
        run_tests, "COVERAGE_HTML", coverage_dir / "html", raising=False
    )

    base_command = [
        sys.executable,
        "-m",
        "pytest",
        "--json-report",
        "--json-report-file",
        str(raw_report_file),
        "--durations=0",
    ]
    monkeypatch.setattr(run_tests, "_BASE_PYTEST_COMMAND", base_command, raising=False)
    monkeypatch.setattr(
        run_tests, "_DEFAULT_TEST_TARGETS", ("tests/unit",), raising=False
    )

    captured_commands: list[list[str]] = []

    def _fake_run_pytest(command: Sequence[str]) -> int:
        captured_commands.append(list(command))
        return 0

    monkeypatch.setattr(run_tests, "run_pytest", _fake_run_pytest)

    structured_template = {
        "meta": {
            "repo": "demo/repo",
            "commit": "abc123",
            "branch": "feature",
            "ts_utc": "2025-01-01T00:00:00+00:00",
            "duration_sec": 12.34,
            "python": "3.11.0",
            "pytest": "8.4.0",
            "exit_code": 0,
        },
        "summary": {
            "total": 20,
            "passed": 19,
            "failed": 1,
            "skipped": 0,
            "xfailed": 0,
            "xpassed": 0,
            "error": 0,
            "success_rate": 0.94,
        },
        "tests": [],
    }

    def _fake_build_structured_report(
        raw: dict[str, Any], exit_code: int
    ) -> dict[str, Any]:
        assert exit_code == 0
        return json.loads(json.dumps(structured_template))

    monkeypatch.setattr(
        run_tests, "build_structured_report", _fake_build_structured_report
    )
    monkeypatch.setattr(run_tests, "_load_raw_report", lambda: {"tests": []})

    captured_configs: list[LoggerConfig] = []

    def _fake_configure_logger(cfg: LoggerConfig) -> object:
        captured_configs.append(cfg)
        return object()

    log_path = tmp_path / "logs" / "run_tests_20240131.log"

    @contextmanager
    def _fake_setup(script_name: str, log_cfg: LoggerConfig, date: str | None = None):
        assert script_name == "run_tests"
        cloned_cfg = LoggerConfig(
            level=log_cfg.level,
            run_id=log_cfg.run_id,
            redact_secrets=log_cfg.redact_secrets,
            stream=log_cfg.stream,
            handlers=list(log_cfg.handlers),
            logger_name=log_cfg.logger_name,
        )
        yield SimpleNamespace(
            log_path=log_path, log_cfg=cloned_cfg, console_stream=None
        )

    monkeypatch.setattr(run_tests, "configure_logger", _fake_configure_logger)
    monkeypatch.setattr(run_tests, "setup_cli_logging", _fake_setup)

    caplog.set_level(logging.ERROR)

    exit_code = run_tests.main(["--date", "20240131"])

    assert exit_code == 1
    assert "below the required 95.00% threshold" in caplog.text
    assert captured_commands, "run_pytest should be invoked"

    payload = json.loads(report_file.read_text(encoding="utf-8"))
    assert payload["summary"]["success_rate"] == pytest.approx(0.94)
    assert captured_configs and captured_configs[-1].level == "INFO"


@pytest.mark.unit
def test_main__writes_reports_when_pytest_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reports_dir = tmp_path / "reports"
    coverage_dir = reports_dir / "coverage"
    raw_report_file = reports_dir / "pytest_raw_report.json"
    report_file = reports_dir / "test_report.json"
    summary_file = reports_dir / "test_summary.md"

    monkeypatch.setattr(run_tests, "ROOT_DIR", tmp_path, raising=False)
    monkeypatch.setattr(run_tests, "REPORTS_DIR", reports_dir, raising=False)
    monkeypatch.setattr(run_tests, "RAW_REPORT_FILE", raw_report_file, raising=False)
    monkeypatch.setattr(run_tests, "REPORT_FILE", report_file, raising=False)
    monkeypatch.setattr(run_tests, "SUMMARY_FILE", summary_file, raising=False)
    monkeypatch.setattr(run_tests, "COVERAGE_DIR", coverage_dir, raising=False)
    monkeypatch.setattr(
        run_tests, "COVERAGE_XML", coverage_dir / "coverage.xml", raising=False
    )
    monkeypatch.setattr(
        run_tests, "COVERAGE_HTML", coverage_dir / "html", raising=False
    )

    base_command = [
        sys.executable,
        "-m",
        "pytest",
        "--json-report",
        "--json-report-file",
        str(raw_report_file),
    ]
    monkeypatch.setattr(run_tests, "_BASE_PYTEST_COMMAND", base_command, raising=False)
    monkeypatch.setattr(
        run_tests, "_DEFAULT_TEST_TARGETS", ("tests/unit",), raising=False
    )
    monkeypatch.setattr(run_tests, "_git_output", lambda *args, **kwargs: "stub")

    captured_commands: list[list[str]] = []

    def _fake_run_pytest(command: Sequence[str]) -> int:
        captured_commands.append(list(command))
        return 2

    monkeypatch.setattr(run_tests, "run_pytest", _fake_run_pytest)
    monkeypatch.setattr(run_tests, "_load_raw_report", lambda: {})

    def _fake_configure_logger(cfg: LoggerConfig) -> object:
        return object()

    monkeypatch.setattr(run_tests, "configure_logger", _fake_configure_logger)

    @contextmanager
    def _fake_setup(script_name: str, log_cfg: LoggerConfig, date: str | None = None):
        yield SimpleNamespace(
            log_path=tmp_path / "logs" / "run_tests.log",
            log_cfg=log_cfg,
            console_stream=None,
        )

    monkeypatch.setattr(run_tests, "setup_cli_logging", _fake_setup)

    exit_code = run_tests.main([])

    assert exit_code == 2
    assert captured_commands, "run_pytest should be invoked"
    assert report_file.exists()
    assert summary_file.exists()

    payload = json.loads(report_file.read_text(encoding="utf-8"))
    summary = payload["summary"]
    assert summary["success_rate"] == 0.0
    assert summary["total"] == 1
    assert summary["error"] == 1

    placeholder = payload["tests"][0]
    assert placeholder["status"] == "error"
    assert "code 2" in placeholder["error"]


@pytest.mark.unit
def test_main__fails_fast_on_summary_filesystem_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reports_dir = tmp_path / "reports"
    coverage_dir = reports_dir / "coverage"
    raw_report_file = reports_dir / "pytest_raw_report.json"
    report_file = reports_dir / "test_report.json"
    summary_file = reports_dir / "test_summary.md"

    monkeypatch.setattr(run_tests, "ROOT_DIR", tmp_path, raising=False)
    monkeypatch.setattr(run_tests, "REPORTS_DIR", reports_dir, raising=False)
    monkeypatch.setattr(run_tests, "RAW_REPORT_FILE", raw_report_file, raising=False)
    monkeypatch.setattr(run_tests, "REPORT_FILE", report_file, raising=False)
    monkeypatch.setattr(run_tests, "SUMMARY_FILE", summary_file, raising=False)
    monkeypatch.setattr(run_tests, "COVERAGE_DIR", coverage_dir, raising=False)
    monkeypatch.setattr(
        run_tests, "COVERAGE_XML", coverage_dir / "coverage.xml", raising=False
    )
    monkeypatch.setattr(
        run_tests, "COVERAGE_HTML", coverage_dir / "html", raising=False
    )

    base_command = [
        sys.executable,
        "-m",
        "pytest",
        "--json-report",
        "--json-report-file",
        str(raw_report_file),
    ]
    monkeypatch.setattr(run_tests, "_BASE_PYTEST_COMMAND", base_command, raising=False)
    monkeypatch.setattr(
        run_tests, "_DEFAULT_TEST_TARGETS", ("tests/unit",), raising=False
    )
    monkeypatch.setattr(run_tests, "_git_output", lambda *args, **kwargs: "stub")

    monkeypatch.setattr(run_tests, "run_pytest", lambda command: 0)
    monkeypatch.setattr(
        run_tests, "_load_raw_report", lambda: {"tests": [], "duration": 0.0}
    )

    def _fake_configure_logger(cfg: LoggerConfig) -> object:
        return object()

    monkeypatch.setattr(run_tests, "configure_logger", _fake_configure_logger)

    log_path = tmp_path / "logs" / "run_tests.log"

    @contextmanager
    def _fake_setup(script_name: str, log_cfg: LoggerConfig, date: str | None = None):
        assert script_name == "run_tests"
        yield SimpleNamespace(log_path=log_path, log_cfg=log_cfg, console_stream=None)

    monkeypatch.setattr(run_tests, "setup_cli_logging", _fake_setup)

    summary_attempts = 0

    def _failing_write_summary(report: dict[str, Any], destination: Path) -> None:
        nonlocal summary_attempts
        summary_attempts += 1
        raise OSError("simulated write failure")

    monkeypatch.setattr(run_tests, "write_summary", _failing_write_summary)

    exit_code = run_tests.main([])

    assert exit_code == run_tests.VALIDATION_FAILURE_EXIT_CODE
    assert summary_attempts == 1
    assert not summary_file.exists()
    assert report_file.exists(), "JSON report should still be produced"
