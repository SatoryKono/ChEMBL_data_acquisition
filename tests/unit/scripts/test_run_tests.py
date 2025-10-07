"""Unit tests for :mod:`scripts.run_tests`."""

from __future__ import annotations

import sys
from collections.abc import Sequence
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pytest

from library.cli import LoggerConfig
from scripts import run_tests


@pytest.mark.unit
def test_run_tests__verbose_creates_debug_log(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
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
    monkeypatch.setattr(run_tests, "COVERAGE_XML", coverage_dir / "coverage.xml", raising=False)
    monkeypatch.setattr(run_tests, "COVERAGE_HTML", coverage_dir / "html", raising=False)

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
    monkeypatch.setattr(run_tests, "_DEFAULT_TEST_TARGETS", ("tests/unit",), raising=False)
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

    exit_code = run_tests.main([
        "--verbose",
        "--date",
        "20240131",
        "--",
        "-k",
        "unit",
    ])

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
