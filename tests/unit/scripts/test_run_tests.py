import json
from types import SimpleNamespace

import pytest

from scripts import run_tests


@pytest.fixture
def _patch_reports(monkeypatch, tmp_path):
    reports_dir = tmp_path / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    report_file = reports_dir / "test_report.json"
    log_file = reports_dir / "test_run.log"
    summary_file = reports_dir / "test_summary.md"

    monkeypatch.setattr(run_tests, "ROOT_DIR", tmp_path)
    monkeypatch.setattr(run_tests, "REPORTS_DIR", reports_dir)
    monkeypatch.setattr(run_tests, "REPORT_FILE", report_file)
    monkeypatch.setattr(run_tests, "LOG_FILE", log_file)
    monkeypatch.setattr(run_tests, "SUMMARY_FILE", summary_file)

    pytest_command = (
        run_tests.sys.executable,
        "-m",
        "pytest",
        "--json-report",
        "--json-report-file",
        str(report_file),
        "--durations=0",
        "-vv",
    )
    monkeypatch.setattr(run_tests, "PYTEST_COMMAND", pytest_command)

    return SimpleNamespace(
        reports_dir=reports_dir,
        report_file=report_file,
        log_file=log_file,
        summary_file=summary_file,
        pytest_command=pytest_command,
    )


def _stub_subprocess(monkeypatch, exit_code):
    calls = {}

    def fake_run(command, stdout=None, stderr=None, check=None):
        calls["command"] = command
        return SimpleNamespace(returncode=exit_code)

    monkeypatch.setattr(run_tests.subprocess, "run", fake_run)
    return calls


def test_main_generates_summary_with_failures(monkeypatch, tmp_path, _patch_reports):
    context = _patch_reports
    report_payload = {
        "summary": {
            "total": 3,
            "passed": 1,
            "failed": 1,
            "error": 1,
            "skipped": 0,
            "xfailed": 0,
            "xpassed": 0,
        },
        "duration": 12.345,
        "tests": [
            {
                "nodeid": "tests/unit/example_test.py::test_failure",
                "outcome": "failed",
                "call": {
                    "longrepr": "AssertionError: expected boom",
                },
            }
        ],
    }
    context.report_file.write_text(json.dumps(report_payload), encoding="utf-8")

    stub_exit_code = 3
    calls = _stub_subprocess(monkeypatch, stub_exit_code)

    run_tests.ensure_reports_directory()
    assert context.reports_dir.is_dir()

    result = run_tests.main()

    assert result == stub_exit_code
    assert calls["command"] == context.pytest_command
    assert context.log_file.exists()
    assert context.summary_file.exists()

    summary_text = context.summary_file.read_text(encoding="utf-8")
    assert "| Metric | Value |" in summary_text
    assert "| Failed | 1 |" in summary_text
    assert "| Errors | 1 |" in summary_text
    assert "| Success rate | 33.3% |" in summary_text
    assert "- Duration: 12.35s" in summary_text
    assert "- `tests/unit/example_test.py::test_failure` (failed)" in summary_text
    assert "AssertionError: expected boom" in summary_text


def test_main_handles_missing_or_corrupt_report(monkeypatch, tmp_path, _patch_reports):
    context = _patch_reports
    context.report_file.write_text("{ not valid json", encoding="utf-8")

    stub_exit_code = 0
    calls = _stub_subprocess(monkeypatch, stub_exit_code)

    run_tests.ensure_reports_directory()
    assert context.reports_dir.exists()

    result = run_tests.main()

    assert result == stub_exit_code
    assert calls["command"] == context.pytest_command
    assert context.log_file.exists()
    assert context.summary_file.exists()

    summary_text = context.summary_file.read_text(encoding="utf-8")
    assert "| Metric | Value |" in summary_text
    assert "| Total tests | 0 |" in summary_text
    assert "All tests passed." in summary_text
