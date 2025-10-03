"""Unit tests for :mod:`tools.run_tests`."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tools import run_tests


def _stub_git_output(_: list[str]) -> str:
    return "stub"


def _make_record(nodeid: str, status: str) -> run_tests.TestRecord:
    return run_tests.TestRecord(
        nodeid=nodeid,
        status=status,
        duration_ms=12.3,
        stdout="",
        stderr="",
    )


@pytest.mark.unit
def test_build_json__stores_fractional_success_rate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(run_tests, "_git_output", _stub_git_output)

    collector = run_tests.ReportCollector()
    collector.start_time = 100.0
    collector.end_time = 100.5
    collector.tests = {
        "tests/test_demo.py::test_pass": _make_record(
            "tests/test_demo.py::test_pass", "passed"
        ),
        "tests/test_demo.py::test_fail": _make_record(
            "tests/test_demo.py::test_fail", "failed"
        ),
    }

    report_path = tmp_path / "report.json"
    payload = run_tests._build_json(collector, report_path=report_path)

    assert report_path.exists()
    saved = json.loads(report_path.read_text(encoding="utf-8"))

    for data in (payload, saved):
        assert data["summary"]["success_rate"] == pytest.approx(0.5)

    summary_path = tmp_path / "summary.md"
    run_tests._write_markdown(summary_path, payload)
    summary_text = summary_path.read_text(encoding="utf-8")

    assert "Success rate: 50.00%" in summary_text


@pytest.mark.unit
def test_build_json__uses_full_success_for_empty_suite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(run_tests, "_git_output", _stub_git_output)

    collector = run_tests.ReportCollector()
    collector.start_time = 50.0
    collector.end_time = 50.0

    report_path = tmp_path / "report.json"
    payload = run_tests._build_json(collector, report_path=report_path)

    assert payload["summary"]["success_rate"] == pytest.approx(1.0)

    summary_path = tmp_path / "summary.md"
    run_tests._write_markdown(summary_path, payload)
    summary_text = summary_path.read_text(encoding="utf-8")

    assert "Success rate: 100.00%" in summary_text
