"""Unit tests for the deprecated :mod:`scripts.run_test_suite` wrapper."""

from __future__ import annotations

import os
from collections.abc import Sequence
from pathlib import Path

import pytest
from scripts import run_test_suite


class _RunTestsSpy:
    def __init__(self) -> None:
        self.calls: list[Sequence[str]] = []
        self.exit_code = 0

    def main(self, args: Sequence[str] | None = None) -> int:
        self.calls.append(tuple(args or []))
        return self.exit_code


@pytest.fixture(autouse=True)
def _reset_hashseed(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("PYTHONHASHSEED", raising=False)


@pytest.mark.unit
def test_main__delegates_to_run_tests_with_default_paths(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    spy = _RunTestsSpy()
    spy.exit_code = 7
    monkeypatch.setattr(run_test_suite.run_tests, "main", spy.main)

    exit_code = run_test_suite.main([])

    assert exit_code == 7
    assert spy.calls == [
        (
            "--json",
            str(Path("reports") / "test_report.json"),
            "--markdown",
            str(Path("reports") / "test_summary.md"),
        )
    ]
    stderr = capsys.readouterr().err.lower()
    assert "deprecated" in stderr


@pytest.mark.unit
def test_main__forwards_pytest_args(monkeypatch: pytest.MonkeyPatch) -> None:
    spy = _RunTestsSpy()
    monkeypatch.setattr(run_test_suite.run_tests, "main", spy.main)

    exit_code = run_test_suite.main(
        ["--pytest-args", "--", "-k", "unit", "tests/unit/test_demo.py"]
    )

    assert exit_code == 0
    assert spy.calls == [
        (
            "--json",
            str(Path("reports") / "test_report.json"),
            "--markdown",
            str(Path("reports") / "test_summary.md"),
            "--",
            "-k",
            "unit",
            "tests/unit/test_demo.py",
        )
    ]


@pytest.mark.unit
def test_main__respects_verbose_and_report_dir(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    spy = _RunTestsSpy()
    monkeypatch.setattr(run_test_suite.run_tests, "main", spy.main)

    exit_code = run_test_suite.main(["--verbose", "--report-dir", str(tmp_path)])

    assert exit_code == 0
    assert spy.calls == [
        (
            "--verbose",
            "--json",
            str(tmp_path / "test_report.json"),
            "--markdown",
            str(tmp_path / "test_summary.md"),
        )
    ]


@pytest.mark.unit
def test_main__sets_pythonhashseed_when_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    spy = _RunTestsSpy()
    monkeypatch.setattr(run_test_suite.run_tests, "main", spy.main)

    run_test_suite.main([])

    assert os.environ["PYTHONHASHSEED"] == "0"


@pytest.mark.unit
def test_main__preserves_existing_pythonhashseed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    spy = _RunTestsSpy()
    monkeypatch.setattr(run_test_suite.run_tests, "main", spy.main)
    monkeypatch.setenv("PYTHONHASHSEED", "123")

    run_test_suite.main([])

    assert os.environ["PYTHONHASHSEED"] == "123"


@pytest.mark.unit
def test_main__emits_warning_for_custom_suite(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    spy = _RunTestsSpy()
    monkeypatch.setattr(run_test_suite.run_tests, "main", spy.main)

    run_test_suite.main(["--suite", "nightly"])

    stderr = capsys.readouterr().err
    assert "deprecated" in stderr.lower()
    assert "nightly" in stderr
