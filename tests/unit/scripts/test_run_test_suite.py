"""Unit tests for the scripts.run_test_suite helper."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from scripts import run_test_suite


def _make_result(name: str, status: str, *, duration: float = 0.1, message: str = "") -> run_test_suite.TestResult:
    return run_test_suite.TestResult(
        name=name,
        status=status,
        duration=duration,
        message=message,
        log_path=None,
    )


@pytest.mark.unit
def test_write_reports__includes_success_rate(tmp_path: Path) -> None:
    plugin = run_test_suite.JsonReportPlugin(tmp_path / "suite.log")
    plugin._results = {
        "tests::passed": _make_result("tests::passed", "passed"),
        "tests::failed": _make_result("tests::failed", "failed", message="boom"),
    }

    results = plugin.results
    summary = run_test_suite.summarize_results(results)

    json_path = tmp_path / "test_report.json"
    run_test_suite._write_json(json_path, results, exit_code=0, summary=summary)
    payload = json.loads(json_path.read_text(encoding="utf-8"))

    assert "success_rate" in payload["summary"]
    assert pytest.approx(payload["summary"]["success_rate"]) == 0.5

    md_path = tmp_path / "test_summary.md"
    run_test_suite._write_summary(md_path, results, exit_code=0, summary=summary)
    summary_text = md_path.read_text(encoding="utf-8")

    assert "* Success rate: 50.00%" in summary_text

    empty_summary = run_test_suite.summarize_results([])
    assert empty_summary["success_rate"] == pytest.approx(1.0)


@pytest.mark.unit
def test_main__returns_error_when_success_rate_below_threshold(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    def fake_pytest_main(pytest_args, plugins):  # type: ignore[override]
        plugin: run_test_suite.JsonReportPlugin = plugins[0]
        log_path = str(plugin._log_path)
        plugin._results = {
            "tests::pass": run_test_suite.TestResult(
                name="tests::pass",
                status="passed",
                duration=0.1,
                message="",
                log_path=log_path,
            ),
            "tests::fail": run_test_suite.TestResult(
                name="tests::fail",
                status="failed",
                duration=0.1,
                message="boom",
                log_path=log_path,
            ),
            "tests::skip": run_test_suite.TestResult(
                name="tests::skip",
                status="skipped",
                duration=0.1,
                message="not run",
                log_path=log_path,
            ),
        }
        return 0

    monkeypatch.setattr(run_test_suite.pytest, "main", fake_pytest_main)

    caplog.set_level(logging.ERROR)
    exit_code = run_test_suite.main(["--report-dir", str(tmp_path), "--suite", "demo"])

    assert run_test_suite.SUCCESS_RATE_THRESHOLD == pytest.approx(0.80)
    assert exit_code == 1
    assert "below the required threshold" in caplog.text

    report_payload = json.loads((tmp_path / "test_report.json").read_text(encoding="utf-8"))
    assert report_payload["summary"]["success_rate"] < run_test_suite.SUCCESS_RATE_THRESHOLD


@pytest.mark.unit
def test_main__passes_when_success_rate_meets_threshold(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    def fake_pytest_main(pytest_args, plugins):  # type: ignore[override]
        plugin: run_test_suite.JsonReportPlugin = plugins[0]
        log_path = str(plugin._log_path)
        plugin._results = {
            "tests::pass": run_test_suite.TestResult(
                name="tests::pass",
                status="passed",
                duration=0.1,
                message="",
                log_path=log_path,
            ),
            "tests::also_pass": run_test_suite.TestResult(
                name="tests::also_pass",
                status="passed",
                duration=0.1,
                message="",
                log_path=log_path,
            ),
            "tests::third_pass": run_test_suite.TestResult(
                name="tests::third_pass",
                status="passed",
                duration=0.1,
                message="",
                log_path=log_path,
            ),
            "tests::fourth_pass": run_test_suite.TestResult(
                name="tests::fourth_pass",
                status="passed",
                duration=0.1,
                message="",
                log_path=log_path,
            ),
            "tests::skip": run_test_suite.TestResult(
                name="tests::skip",
                status="skipped",
                duration=0.1,
                message="not run",
                log_path=log_path,
            ),
        }
        return 0

    monkeypatch.setattr(run_test_suite.pytest, "main", fake_pytest_main)

    caplog.set_level(logging.ERROR)
    exit_code = run_test_suite.main(["--report-dir", str(tmp_path), "--suite", "demo"])

    assert exit_code == 0
    assert "below the required threshold" not in caplog.text

    report_payload = json.loads((tmp_path / "test_report.json").read_text(encoding="utf-8"))
    assert report_payload["summary"]["success_rate"] >= run_test_suite.SUCCESS_RATE_THRESHOLD
