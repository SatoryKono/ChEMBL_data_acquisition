"""Unit tests for the scripts.run_test_suite helper."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Callable

import pytest

from scripts import run_test_suite


def _make_result(
    name: str,
    status: str,
    *,
    duration_ms: float = 100.0,
    stdout: str = "",
    stderr: str = "",
    error: str | None = None,
    log: list[str] | None = None,
) -> run_test_suite.TestResult:
    return run_test_suite.TestResult(
        nodeid=name,
        status=status,
        duration_ms=duration_ms,
        stdout=stdout,
        stderr=stderr,
        log=list(log or []),
        error=error,
    )


def _patch_meta(monkeypatch: pytest.MonkeyPatch, duration_factory: Callable[[float], dict[str, object]]) -> None:
    monkeypatch.setattr(run_test_suite, "_gather_meta", duration_factory)


@pytest.mark.unit
def test_write_reports__includes_success_rate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _patch_meta(
        monkeypatch,
        lambda duration: {
            "repo": "demo/repo",
            "commit": "deadbeef",
            "branch": "main",
            "ts_utc": "2020-01-01T00:00:00Z",
            "duration_sec": round(duration, 3),
            "python": "3.11.0",
            "pytest": "8.0.0",
        },
    )

    results = [
        _make_result("tests::passed", "passed"),
        _make_result("tests::failed", "failed", error="boom"),
    ]

    payload = run_test_suite._build_report(results, duration_sec=12.345)

    json_path = tmp_path / "test_report.json"
    run_test_suite._write_json(json_path, payload)
    saved = json.loads(json_path.read_text(encoding="utf-8"))

    assert saved["summary"]["success_rate"] == pytest.approx(0.5)
    assert saved["meta"]["repo"] == "demo/repo"

    md_path = tmp_path / "test_summary.md"
    run_test_suite._write_summary(md_path, payload)
    summary_text = md_path.read_text(encoding="utf-8")

    assert "- Success rate: 50.00%" in summary_text
    assert "- `tests::failed`: `boom`" in summary_text

    empty_summary = run_test_suite.summarize_results([])
    assert empty_summary["success_rate"] == pytest.approx(1.0)


@pytest.mark.unit
def test_main__returns_error_when_success_rate_below_threshold(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    _patch_meta(
        monkeypatch,
        lambda duration: {
            "repo": "demo/repo",
            "commit": "deadbeef",
            "branch": "main",
            "ts_utc": "2020-01-01T00:00:00Z",
            "duration_sec": round(duration, 3),
            "python": "3.11.0",
            "pytest": "8.0.0",
        },
    )

    def fake_pytest_main(pytest_args, plugins):  # type: ignore[override]
        plugin: run_test_suite.JsonReportPlugin = plugins[0]
        plugin._started_at = 0.0
        plugin._finished_at = 1.0
        plugin._results = {
            **{
                f"tests::pass_{index}": _make_result(f"tests::pass_{index}", "passed")
                for index in range(18)
            },
            "tests::fail": _make_result("tests::fail", "failed", error="boom"),
        }
        return 0

    monkeypatch.setattr(run_test_suite.pytest, "main", fake_pytest_main)

    caplog.set_level(logging.ERROR)
    exit_code = run_test_suite.main(["--report-dir", str(tmp_path), "--suite", "demo"])

    assert run_test_suite.SUCCESS_RATE_THRESHOLD == pytest.approx(0.95)
    assert exit_code == 1
    assert "below the required threshold" in caplog.text

    report_payload = json.loads((tmp_path / "test_report.json").read_text(encoding="utf-8"))
    assert report_payload["summary"]["success_rate"] < run_test_suite.SUCCESS_RATE_THRESHOLD


@pytest.mark.unit
def test_main__passes_when_success_rate_meets_threshold(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    _patch_meta(
        monkeypatch,
        lambda duration: {
            "repo": "demo/repo",
            "commit": "deadbeef",
            "branch": "main",
            "ts_utc": "2020-01-01T00:00:00Z",
            "duration_sec": round(duration, 3),
            "python": "3.11.0",
            "pytest": "8.0.0",
        },
    )

    def fake_pytest_main(pytest_args, plugins):  # type: ignore[override]
        plugin: run_test_suite.JsonReportPlugin = plugins[0]
        plugin._started_at = 0.0
        plugin._finished_at = 1.0
        plugin._results = {
            **{
                f"tests::pass_{index}": _make_result(f"tests::pass_{index}", "passed")
                for index in range(19)
            },
            "tests::skip": _make_result("tests::skip", "skipped", error="not run"),
        }
        return 0

    monkeypatch.setattr(run_test_suite.pytest, "main", fake_pytest_main)

    caplog.set_level(logging.ERROR)
    exit_code = run_test_suite.main(["--report-dir", str(tmp_path), "--suite", "demo"])

    assert exit_code == 0
    assert "below the required threshold" not in caplog.text

    report_payload = json.loads((tmp_path / "test_report.json").read_text(encoding="utf-8"))
    assert report_payload["summary"]["success_rate"] >= run_test_suite.SUCCESS_RATE_THRESHOLD
