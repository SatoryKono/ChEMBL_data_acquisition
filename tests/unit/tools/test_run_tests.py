"""Unit tests for :mod:`tests.run_tests`."""

from __future__ import annotations

import importlib.util
import io
import json
import sys
from collections.abc import Sequence
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any, Protocol, TYPE_CHECKING

import pytest


_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_RUN_TESTS_PATH = _PROJECT_ROOT / "tests" / "run_tests.py"


if TYPE_CHECKING:
    from tests.run_tests import LoggerConfig as _LoggerConfig
    from tests.run_tests import ReportCollector as _ReportCollector
    from tests.run_tests import TestRecord as _TestRecord
else:  # pragma: no cover - used only for typing
    _LoggerConfig = Any  # type: ignore[assignment]
    _ReportCollector = Any  # type: ignore[assignment]
    _TestRecord = Any  # type: ignore[assignment]


class RunTestsModule(Protocol):
    """Subset of :mod:`tests.run_tests` used within this test module."""

    LoggerConfig: type[_LoggerConfig]
    TestRecord: type[_TestRecord]
    ReportCollector: type[_ReportCollector]
    pytest: ModuleType
    DEPRECATION_MESSAGE: str

    def _build_json(
        self, collector: _ReportCollector, *, report_path: Path
    ) -> dict[str, Any]:
        ...

    def _write_markdown(self, path: Path, payload: dict[str, Any]) -> None:
        ...

    def main(self, argv: Sequence[str]) -> int:
        ...


LoggerConfig = _LoggerConfig
ReportCollector = _ReportCollector
TestRecord = _TestRecord


@pytest.fixture(name="run_tests_module")
def fixture_run_tests_module() -> RunTestsModule:
    """Load :mod:`tests.run_tests` in isolation and restore ``sys.modules`` afterwards."""

    module_name = "tests.run_tests"
    spec = importlib.util.spec_from_file_location(module_name, _RUN_TESTS_PATH)
    if spec is None or spec.loader is None:  # pragma: no cover - defensive guard
        raise RuntimeError(f"Unable to load run_tests module from {_RUN_TESTS_PATH}")

    module = importlib.util.module_from_spec(spec)
    previous = sys.modules.get(module_name)
    sys.modules[module_name] = module

    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(module_name, None)
        if previous is not None:
            sys.modules[module_name] = previous
        raise

    try:
        yield module
    finally:
        sys.modules.pop(module_name, None)
        if previous is not None:
            sys.modules[module_name] = previous


def _stub_git_output(_: list[str]) -> str:
    return "stub"


def _make_record(run_tests: RunTestsModule, nodeid: str, status: str) -> TestRecord:
    return run_tests.TestRecord(
        nodeid=nodeid,
        status=status,
        duration_ms=12.3,
        stdout="",
        stderr="",
    )


@pytest.mark.unit
def test_build_json__stores_fractional_success_rate(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    run_tests_module: RunTestsModule,
) -> None:
    _install_fake_cli_logging(run_tests_module, monkeypatch, tmp_path)
    monkeypatch.setattr(run_tests_module, "_git_output", _stub_git_output)

    collector = run_tests_module.ReportCollector()
    collector.start_time = 100.0
    collector.end_time = 100.5
    collector.tests = {
        "tests/test_demo.py::test_pass": _make_record(
            run_tests_module,
            "tests/test_demo.py::test_pass", "passed"
        ),
        "tests/test_demo.py::test_fail": _make_record(
            run_tests_module,
            "tests/test_demo.py::test_fail", "failed"
        ),
    }

    report_path = tmp_path / "report.json"
    payload = run_tests_module._build_json(collector, report_path=report_path)

    assert report_path.exists()
    saved = json.loads(report_path.read_text(encoding="utf-8"))

    for data in (payload, saved):
        assert data["summary"]["success_rate"] == pytest.approx(0.5)

    summary_path = tmp_path / "summary.md"
    run_tests_module._write_markdown(summary_path, payload)
    summary_text = summary_path.read_text(encoding="utf-8")

    assert "Success rate: 50.00%" in summary_text


@pytest.mark.unit
def test_build_json__counts_xfailed_and_excludes_skipped(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    run_tests_module: RunTestsModule,
) -> None:
    _install_fake_cli_logging(run_tests_module, monkeypatch, tmp_path)
    monkeypatch.setattr(run_tests_module, "_git_output", _stub_git_output)

    collector = run_tests_module.ReportCollector()
    collector.start_time = 10.0
    collector.end_time = 10.0
    collector.tests = {
        "tests/test_demo.py::test_pass": _make_record(
            run_tests_module, "tests/test_demo.py::test_pass", "passed"
        ),
        "tests/test_demo.py::test_skip": _make_record(
            run_tests_module, "tests/test_demo.py::test_skip", "skipped"
        ),
        "tests/test_demo.py::test_xfail": _make_record(
            run_tests_module, "tests/test_demo.py::test_xfail", "xfailed"
        ),
        "tests/test_demo.py::test_fail": _make_record(
            run_tests_module, "tests/test_demo.py::test_fail", "failed"
        ),
        "tests/test_demo.py::test_error": _make_record(
            run_tests_module, "tests/test_demo.py::test_error", "error"
        ),
        "tests/test_demo.py::test_xpass": _make_record(
            run_tests_module, "tests/test_demo.py::test_xpass", "xpassed"
        ),
    }
    collector.tests["tests/test_demo.py::test_error"].error = "boom"

    report_path = tmp_path / "report.json"
    payload = run_tests_module._build_json(collector, report_path=report_path)

    summary = payload["summary"]
    assert summary["total"] == 6
    assert summary["passed"] == 1
    assert summary["xfailed"] == 1
    assert summary["skipped"] == 1
    assert summary["failed"] == 1
    assert summary["xpassed"] == 1
    assert summary["error"] == 1
    assert summary["success_rate"] == pytest.approx(0.4)

    summary_path = tmp_path / "summary.md"
    run_tests_module._write_markdown(summary_path, payload)
    summary_text = summary_path.read_text(encoding="utf-8")
    assert "Success rate: 40.00%" in summary_text


@pytest.mark.unit
def test_build_json__uses_full_success_for_empty_suite(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    run_tests_module: RunTestsModule,
) -> None:
    _install_fake_cli_logging(run_tests_module, monkeypatch, tmp_path)
    monkeypatch.setattr(run_tests_module, "_git_output", _stub_git_output)

    collector = run_tests_module.ReportCollector()
    collector.start_time = 50.0
    collector.end_time = 50.0

    report_path = tmp_path / "report.json"
    payload = run_tests_module._build_json(collector, report_path=report_path)

    assert payload["summary"]["success_rate"] == pytest.approx(1.0)

    summary_path = tmp_path / "summary.md"
    run_tests_module._write_markdown(summary_path, payload)
    summary_text = summary_path.read_text(encoding="utf-8")

    assert "Success rate: 100.00%" in summary_text


@pytest.mark.unit
def test_main__downgrades_exit_code_when_threshold_not_met(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    run_tests_module: RunTestsModule,
) -> None:
    log_path, captured_cfgs = _install_fake_cli_logging(
        run_tests_module, monkeypatch, tmp_path
    )
    monkeypatch.setattr(run_tests_module, "_git_output", _stub_git_output)

    def _fake_pytest_main(pytest_args: list[str], plugins: list[object]) -> int:
        assert plugins, "Collector plugin must be provided"
        collector = plugins[0]
        assert isinstance(collector, run_tests_module.ReportCollector)
        assert captured_cfgs, "configure_logger should be called before pytest executes"
        assert f"--log-file={log_path}" in pytest_args
        assert f"--log-file-level={captured_cfgs[-1].level}" in pytest_args

        collector.start_time = 0.0
        collector.end_time = 0.5
        collector.tests = {
            "tests/test_demo.py::test_pass": _make_record(
                run_tests_module,
                "tests/test_demo.py::test_pass", "passed"
            ),
            "tests/test_demo.py::test_fail": _make_record(
                run_tests_module,
                "tests/test_demo.py::test_fail", "failed"
            ),
        }
        collector.tests["tests/test_demo.py::test_fail"].error = "boom"
        return 0

    monkeypatch.setattr(run_tests_module.pytest, "main", _fake_pytest_main)

    json_path = tmp_path / "report.json"
    markdown_path = tmp_path / "summary.md"

    exit_code = run_tests_module.main([
        "--json",
        str(json_path),
        "--markdown",
        str(markdown_path),
    ])

    assert exit_code == 1

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["summary"]["success_rate"] == pytest.approx(0.5)

    markdown_text = markdown_path.read_text(encoding="utf-8")
    assert "Success rate: 50.00%" in markdown_text
    assert captured_cfgs[-1].level == "INFO"


@pytest.mark.unit
def test_main__verbose_propagates_debug_logging(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    run_tests_module: RunTestsModule,
) -> None:
    log_path, captured_cfgs = _install_fake_cli_logging(
        run_tests_module, monkeypatch, tmp_path
    )
    monkeypatch.setattr(run_tests_module, "_git_output", _stub_git_output)

    def _fake_pytest_main(pytest_args: list[str], plugins: list[object]) -> int:
        assert f"--log-file={log_path}" in pytest_args
        assert f"--log-file-level=DEBUG" in pytest_args
        return 0

    monkeypatch.setattr(run_tests_module.pytest, "main", _fake_pytest_main)

    exit_code = run_tests_module.main([
        "--json",
        str(tmp_path / "report.json"),
        "--markdown",
        str(tmp_path / "summary.md"),
        "--verbose",
    ])

    assert exit_code == 0
    assert captured_cfgs and captured_cfgs[-1].level == "DEBUG"


@pytest.mark.unit
def test_main__emits_deprecation_warning(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    run_tests_module: RunTestsModule,
) -> None:
    log_path, _ = _install_fake_cli_logging(run_tests_module, monkeypatch, tmp_path)
    monkeypatch.setattr(run_tests_module, "_git_output", _stub_git_output)

    def _fake_pytest_main(pytest_args: list[str], plugins: list[object]) -> int:
        assert f"--log-file={log_path}" in pytest_args
        return 0

    monkeypatch.setattr(run_tests_module.pytest, "main", _fake_pytest_main)

    with pytest.warns(DeprecationWarning) as recorded:
        exit_code = run_tests_module.main(
            [
                "--json",
                str(tmp_path / "report.json"),
                "--markdown",
                str(tmp_path / "summary.md"),
            ]
        )

    assert exit_code == 0
    assert recorded
    assert run_tests_module.DEPRECATION_MESSAGE in str(recorded[0].message)


def _install_fake_cli_logging(
    run_tests: RunTestsModule, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> tuple[Path, list[LoggerConfig]]:
    log_path = tmp_path / "run_tests.log"
    captured_cfgs: list[LoggerConfig] = []

    def _fake_configure(cfg: LoggerConfig, *_: object, **__: object) -> object:
        captured_cfgs.append(cfg)
        return object()

    @contextmanager
    def _fake_setup(script_name: str, log_cfg: LoggerConfig, *_: object, **__: object):
        assert script_name == "run_tests"
        stream = io.StringIO()
        cloned_cfg = run_tests.LoggerConfig(
            level=log_cfg.level,
            run_id=log_cfg.run_id,
            redact_secrets=log_cfg.redact_secrets,
            stream=stream,
            handlers=list(log_cfg.handlers),
            logger_name=log_cfg.logger_name,
        )
        yield SimpleNamespace(log_path=log_path, log_cfg=cloned_cfg, console_stream=stream)

    monkeypatch.setattr(run_tests, "configure_logger", _fake_configure)
    monkeypatch.setattr(run_tests, "setup_cli_logging", _fake_setup)
    return log_path, captured_cfgs
