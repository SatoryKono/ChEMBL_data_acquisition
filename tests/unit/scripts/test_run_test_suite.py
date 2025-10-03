"""Tests for :mod:`scripts.run_test_suite`."""

from __future__ import annotations

import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import run_test_suite


def _make_report(
    nodeid: str,
    outcome: str,
    duration: float,
    longrepr: object = "",
) -> SimpleNamespace:
    """Create a lightweight object mimicking :class:`pytest.TestReport`."""

    return SimpleNamespace(
        nodeid=nodeid,
        outcome=outcome,
        duration=duration,
        longrepr=longrepr,
    )


def test_run_test_suite_creates_reports(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure the CLI writes JSON/Markdown reports and configures pytest correctly."""

    plugin = run_test_suite.JsonReportPlugin(Path("placeholder.log"))

    def fake_plugin_factory(log_path: Path) -> run_test_suite.JsonReportPlugin:
        plugin._log_path = log_path  # type: ignore[attr-defined]
        return plugin

    captured: dict[str, object] = {}

    def fake_pytest_main(args: list[str], plugins: list[run_test_suite.JsonReportPlugin]) -> int:
        captured["args"] = args
        captured["plugins"] = plugins
        reports = [
            _make_report("tests/test_module.py::test_pass", "passed", 0.75),
            _make_report(
                "tests/test_module.py::test_fail",
                "failed",
                1.50,
                longrepr="AssertionError: boom",
            ),
            _make_report(
                "tests/test_module.py::test_skip",
                "skipped",
                0.25,
                longrepr=("tests/test_module.py", 42, "because reasons"),
            ),
        ]
        for report in reports:
            plugins[0].pytest_runtest_logreport(report)
        return 0

    monkeypatch.setattr(run_test_suite, "JsonReportPlugin", fake_plugin_factory)
    monkeypatch.setattr(run_test_suite.pytest, "main", fake_pytest_main)
    monkeypatch.delenv("PYTHONHASHSEED", raising=False)

    exit_code = run_test_suite.main([
        "--suite",
        "unit",
        "--report-dir",
        str(tmp_path),
    ])

    assert exit_code == 0

    expected_log_path = tmp_path / "logs" / "unit.log"
    assert captured["args"] == [
        "tests",
        f"--log-file={expected_log_path}",
        "--log-file-level=DEBUG",
        "--maxfail=0",
    ]
    assert captured["plugins"] == [plugin]
    assert os.environ["PYTHONHASHSEED"] == "0"

    report_path = tmp_path / "test_report.json"
    summary_path = tmp_path / "test_summary.md"
    assert report_path.exists()
    assert summary_path.exists()

    with report_path.open(encoding="utf-8") as fh:
        report_payload = json.load(fh)

    assert report_payload["exit_code"] == 0
    tests_by_name = {entry["name"]: entry for entry in report_payload["tests"]}
    assert set(tests_by_name) == {
        "tests/test_module.py::test_pass",
        "tests/test_module.py::test_fail",
        "tests/test_module.py::test_skip",
    }

    assert tests_by_name["tests/test_module.py::test_pass"]["status"] == "passed"
    assert tests_by_name["tests/test_module.py::test_pass"]["duration"] == pytest.approx(0.75)
    assert tests_by_name["tests/test_module.py::test_pass"]["message"] == ""
    assert tests_by_name["tests/test_module.py::test_pass"]["log_path"] == str(expected_log_path)

    assert tests_by_name["tests/test_module.py::test_fail"]["status"] == "failed"
    assert tests_by_name["tests/test_module.py::test_fail"]["duration"] == pytest.approx(1.50)
    assert tests_by_name["tests/test_module.py::test_fail"]["message"] == "AssertionError: boom"
    assert tests_by_name["tests/test_module.py::test_fail"]["log_path"] == str(expected_log_path)

    assert tests_by_name["tests/test_module.py::test_skip"]["status"] == "skipped"
    assert tests_by_name["tests/test_module.py::test_skip"]["duration"] == pytest.approx(0.25)
    assert tests_by_name["tests/test_module.py::test_skip"]["message"] == "because reasons"
    assert tests_by_name["tests/test_module.py::test_skip"]["log_path"] == str(expected_log_path)

    summary_text = summary_path.read_text(encoding="utf-8")
    assert "# Test suite summary" in summary_text
    assert "* Exit code: 0" in summary_text
    assert "* Tests collected: 3" in summary_text
    assert "## Failing tests" in summary_text
    assert "- `tests/test_module.py::test_fail` (1.50s)" in summary_text
    assert "AssertionError: boom" in summary_text
    assert "## Skipped tests" in summary_text
    assert "- `tests/test_module.py::test_skip` (0.25s)" in summary_text
    assert "Reason: because reasons" in summary_text


def test_run_test_suite_propagates_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Non-zero pytest exit codes should be forwarded and reported as failures."""

    plugin = run_test_suite.JsonReportPlugin(Path("placeholder.log"))

    def fake_plugin_factory(log_path: Path) -> run_test_suite.JsonReportPlugin:
        plugin._log_path = log_path  # type: ignore[attr-defined]
        return plugin

    failing_report = _make_report(
        "tests/test_module.py::test_fail",
        "failed",
        2.0,
        longrepr="RuntimeError: busted",
    )

    def fake_pytest_main(args: list[str], plugins: list[run_test_suite.JsonReportPlugin]) -> int:
        plugins[0].pytest_runtest_logreport(failing_report)
        return 2

    monkeypatch.setattr(run_test_suite, "JsonReportPlugin", fake_plugin_factory)
    monkeypatch.setattr(run_test_suite.pytest, "main", fake_pytest_main)
    monkeypatch.delenv("PYTHONHASHSEED", raising=False)

    exit_code = run_test_suite.main([
        "--suite",
        "unit",
        "--report-dir",
        str(tmp_path),
    ])

    assert exit_code == 2

    summary_text = (tmp_path / "test_summary.md").read_text(encoding="utf-8")
    assert "* Exit code: 2" in summary_text
    assert "- `tests/test_module.py::test_fail` (2.00s)" in summary_text

    report_payload = json.loads((tmp_path / "test_report.json").read_text(encoding="utf-8"))
    assert report_payload["exit_code"] == 2
    tests_by_name = {entry["name"]: entry for entry in report_payload["tests"]}
    assert tests_by_name["tests/test_module.py::test_fail"]["status"] == "failed"

