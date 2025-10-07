"""Unit tests for the scripts.run_test_suite helper."""

from __future__ import annotations

import io
import json
import logging
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

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
    log_path = tmp_path / "suite.log"

    captured_cfgs: list[run_test_suite.LoggerConfig] = []

    def fake_configure_logger(cfg: run_test_suite.LoggerConfig, *_: object, **__: object) -> object:
        captured_cfgs.append(cfg)
        return object()

    @contextmanager
    def fake_setup_cli_logging(script_name: str, log_cfg: run_test_suite.LoggerConfig, *_: object, **__: object):
        assert script_name == "run_test_suite"
        stream = io.StringIO()
        cloned_cfg = run_test_suite.LoggerConfig(
            level=log_cfg.level,
            run_id=log_cfg.run_id,
            redact_secrets=log_cfg.redact_secrets,
            stream=stream,
            handlers=list(log_cfg.handlers),
            logger_name=log_cfg.logger_name,
        )
        yield SimpleNamespace(log_path=log_path, log_cfg=cloned_cfg, console_stream=stream)

    monkeypatch.setattr(run_test_suite, "configure_logger", fake_configure_logger)
    monkeypatch.setattr(run_test_suite, "setup_cli_logging", fake_setup_cli_logging)

    def fake_pytest_main(pytest_args, plugins):  # type: ignore[override]
        plugin: run_test_suite.JsonReportPlugin = plugins[0]
        log_path = str(plugin._log_path)
        assert captured_cfgs, "configure_logger should be called before pytest executes"
        assert f"--log-file={log_path}" in pytest_args
        assert f"--log-file-level={captured_cfgs[-1].level}" in pytest_args
        plugin._results = {
            **{
                f"tests::pass_{index}": run_test_suite.TestResult(
                    name=f"tests::pass_{index}",
                    status="passed",
                    duration=0.1,
                    message="",
                    log_path=log_path,
                )
                for index in range(18)
            },
            "tests::fail": run_test_suite.TestResult(
                name="tests::fail",
                status="failed",
                duration=0.1,
                message="boom",
                log_path=log_path,
            ),
        }
        return 0

    monkeypatch.setattr(run_test_suite.pytest, "main", fake_pytest_main)

    caplog.set_level(logging.ERROR)
    exit_code = run_test_suite.main(["--report-dir", str(tmp_path), "--suite", "demo"])

    assert run_test_suite.SUCCESS_RATE_THRESHOLD == pytest.approx(0.95)
    assert exit_code == 1
    assert "below the required threshold" in caplog.text
    assert captured_cfgs and captured_cfgs[-1].level == "INFO"
    assert captured_cfgs[-1].logger_name == "run_test_suite"

    report_payload = json.loads((tmp_path / "test_report.json").read_text(encoding="utf-8"))
    assert report_payload["summary"]["success_rate"] < run_test_suite.SUCCESS_RATE_THRESHOLD


@pytest.mark.unit
def test_main__passes_when_success_rate_meets_threshold(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    log_path = tmp_path / "suite.log"
    captured_cfgs: list[run_test_suite.LoggerConfig] = []

    def fake_configure_logger(cfg: run_test_suite.LoggerConfig, *_: object, **__: object) -> object:
        captured_cfgs.append(cfg)
        return object()

    @contextmanager
    def fake_setup_cli_logging(script_name: str, log_cfg: run_test_suite.LoggerConfig, *_: object, **__: object):
        assert script_name == "run_test_suite"
        stream = io.StringIO()
        cloned_cfg = run_test_suite.LoggerConfig(
            level=log_cfg.level,
            run_id=log_cfg.run_id,
            redact_secrets=log_cfg.redact_secrets,
            stream=stream,
            handlers=list(log_cfg.handlers),
            logger_name=log_cfg.logger_name,
        )
        yield SimpleNamespace(log_path=log_path, log_cfg=cloned_cfg, console_stream=stream)

    monkeypatch.setattr(run_test_suite, "configure_logger", fake_configure_logger)
    monkeypatch.setattr(run_test_suite, "setup_cli_logging", fake_setup_cli_logging)

    def fake_pytest_main(pytest_args, plugins):  # type: ignore[override]
        plugin: run_test_suite.JsonReportPlugin = plugins[0]
        log_path = str(plugin._log_path)
        assert captured_cfgs, "configure_logger should be called before pytest executes"
        assert f"--log-file={log_path}" in pytest_args
        plugin._results = {
            **{
                f"tests::pass_{index}": run_test_suite.TestResult(
                    name=f"tests::pass_{index}",
                    status="passed",
                    duration=0.1,
                    message="",
                    log_path=log_path,
                )
                for index in range(19)
            },
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
    exit_code = run_test_suite.main(["--report-dir", str(tmp_path), "--suite", "demo", "--verbose"])

    assert exit_code == 0
    assert "below the required threshold" not in caplog.text
    assert captured_cfgs and captured_cfgs[-1].level == "DEBUG"

    report_payload = json.loads((tmp_path / "test_report.json").read_text(encoding="utf-8"))
    assert report_payload["summary"]["success_rate"] >= run_test_suite.SUCCESS_RATE_THRESHOLD
