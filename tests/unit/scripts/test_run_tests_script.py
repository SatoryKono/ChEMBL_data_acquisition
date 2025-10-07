"""Unit tests for :mod:`scripts.run_tests`."""

from __future__ import annotations

from contextlib import contextmanager
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Iterable, Sequence

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

    tests_root = tmp_path / "tests"
    tests_root.mkdir()

    monkeypatch.setattr(run_tests, "ROOT_DIR", tmp_path)
    monkeypatch.setattr(run_tests, "REPORTS_DIR", reports_dir, raising=False)
    monkeypatch.setattr(run_tests, "RAW_REPORT_FILE", raw_report_file, raising=False)
    monkeypatch.setattr(run_tests, "REPORT_FILE", report_file, raising=False)
    monkeypatch.setattr(run_tests, "SUMMARY_FILE", summary_file, raising=False)
    monkeypatch.setattr(run_tests, "COVERAGE_DIR", coverage_dir, raising=False)
    monkeypatch.setattr(run_tests, "COVERAGE_XML", coverage_dir / "coverage.xml", raising=False)
    monkeypatch.setattr(run_tests, "COVERAGE_HTML", coverage_dir / "html", raising=False)
    monkeypatch.setattr(run_tests, "TESTS_ROOT", tests_root, raising=False)

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
        run_tests,
        "_discover_test_targets",
        lambda *, exclude=None: ("tests/unit",),
        raising=False,
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
        "success_rate": 0.5,
    }

    failure_entry = next(
        item
        for item in structured["tests"]
        if item["nodeid"].endswith("test_failure")
    )
    assert failure_entry["status"] == "failed"
    assert failure_entry["error"] == "AssertionError: boom\nline 1\nline 2"


@pytest.mark.unit
def test_build_summary_markdown__renders_error_messages_from_json() -> None:
    report = {
        "meta": {
            "repo": "demo/repo",
            "commit": "abc123",
            "branch": "feature",
            "ts_utc": "2025-01-01T00:00:00+00:00",
            "duration_sec": 3.14,
        },
        "summary": {
            "total": 1,
            "passed": 0,
            "failed": 1,
            "skipped": 0,
            "xfailed": 0,
            "xpassed": 0,
            "error": 0,
            "success_rate": 0.0,
        },
        "tests": [
            {
                "nodeid": "tests/unit/test_example.py::test_failure",
                "status": "failed",
                "duration_ms": 123.0,
                "stdout": "",
                "stderr": "",
                "log": [],
                "error": "AssertionError: boom\nline 1\nline 2",
            }
        ],
    }

    summary_md = run_tests.build_summary_markdown(report)

    assert "- `tests/unit/test_example.py::test_failure` (failed)" in summary_md
    assert summary_md.count("```") == 2
    assert "AssertionError: boom" in summary_md
    assert "line 1" in summary_md and "line 2" in summary_md


@pytest.mark.unit
def test_discover_test_targets__captures_new_directories(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tests_root = tmp_path / "tests"
    unit_dir = tests_root / "unit"
    new_suite = tests_root / "new_suite"
    helpers_dir = tests_root / "helpers"

    for directory in (unit_dir, new_suite, helpers_dir):
        directory.mkdir(parents=True)

    (unit_dir / "test_example.py").write_text("def test_unit():\n    assert True\n", encoding="utf-8")
    (new_suite / "test_new.py").write_text("def test_new():\n    assert True\n", encoding="utf-8")
    (helpers_dir / "util.py").write_text("def helper():\n    return 1\n", encoding="utf-8")

    monkeypatch.setattr(run_tests, "TESTS_ROOT", tests_root, raising=False)

    discovered = run_tests._discover_test_targets()

    expected_targets = {
        str(unit_dir),
        str(new_suite),
        str(helpers_dir),
    }
    assert set(discovered) == expected_targets

    excluded = run_tests._discover_test_targets(exclude=["helpers"])
    assert str(helpers_dir) not in excluded
    assert str(new_suite) in excluded
    assert str(unit_dir) in excluded


@pytest.mark.unit
def test_main__forwards_exclude_arguments(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tests_root = tmp_path / "tests"
    tests_root.mkdir()

    reports_dir = tmp_path / "reports"
    coverage_dir = reports_dir / "coverage"
    raw_report_file = reports_dir / "pytest_raw_report.json"
    report_file = reports_dir / "test_report.json"
    summary_file = reports_dir / "test_summary.md"

    monkeypatch.setattr(run_tests, "ROOT_DIR", tmp_path)
    monkeypatch.setattr(run_tests, "TESTS_ROOT", tests_root, raising=False)
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

    captured_exclude: list[str] | None = None

    def _fake_discover(*, exclude: Iterable[str] | None = None) -> tuple[str, ...]:
        nonlocal captured_exclude
        captured_exclude = list(exclude or [])
        return ("tests/unit",)

    monkeypatch.setattr(run_tests, "_discover_test_targets", _fake_discover)
    monkeypatch.setattr(run_tests, "_git_output", lambda *args, **kwargs: "stub")
    monkeypatch.setattr(run_tests, "run_pytest", lambda command: 0)
    monkeypatch.setattr(run_tests, "configure_logger", lambda cfg: object())

    @contextmanager
    def _fake_setup(script_name: str, log_cfg: LoggerConfig, date: str | None = None):
        yield SimpleNamespace(
            log_path=tmp_path / "log.txt",
            log_cfg=log_cfg,
            console_stream=None,
        )

    monkeypatch.setattr(run_tests, "setup_cli_logging", _fake_setup)
    monkeypatch.setattr(
        run_tests,
        "_load_raw_report",
        lambda: {"tests": [], "duration": 0.0},
    )

    exit_code = run_tests.main(["--exclude", "helpers", "--", "-k", "unit"])

    assert exit_code == 0
    assert captured_exclude == ["helpers"]
