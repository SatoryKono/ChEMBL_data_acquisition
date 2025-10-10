"""CLI tests for :mod:`tools.make_md_summary`."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Callable

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.unit
def test_make_md_summary__creates_markdown_from_valid_report(tmp_path: Path) -> None:
    """The CLI should transform a valid JSON report into a Markdown summary."""

    report_path = tmp_path / "report.json"
    output_path = tmp_path / "summary.md"

    report_payload = {
        "meta": {
            "repo": "SatoryKono/ChEMBL_data_acquisition",
            "commit": "abc1234",
            "branch": "feature/test",
            "ts_utc": "2024-01-01T00:00:00+00:00",
            "duration_sec": 12.5,
        },
        "summary": {
            "total": 2,
            "passed": 2,
            "failed": 0,
            "skipped": 0,
            "xfailed": 0,
            "xpassed": 0,
            "error": 0,
            "success_rate": 1.0,
        },
        "tests": [
            {
                "nodeid": "tests/unit/test_example.py::test_one",
                "status": "passed",
                "duration_ms": 5.0,
                "stdout": "",
                "stderr": "",
                "log": [],
                "error": None,
            },
            {
                "nodeid": "tests/unit/test_example.py::test_two",
                "status": "passed",
                "duration_ms": 7.0,
                "stdout": "",
                "stderr": "",
                "log": [],
                "error": None,
            },
        ],
    }
    report_path.write_text(json.dumps(report_payload), encoding="utf-8")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "tools.make_md_summary",
            "--input",
            str(report_path),
            "--output",
            str(output_path),
        ],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert result.stdout == ""
    assert result.stderr == ""
    assert output_path.exists()
    summary_text = output_path.read_text(encoding="utf-8")
    assert "# Test Summary" in summary_text
    assert "Success rate: 100.00%" in summary_text
    assert "|     2 |     2 |     0" in summary_text


@pytest.mark.unit
@pytest.mark.parametrize(
    "input_factory",
    (
        pytest.param(
            lambda directory: _create_invalid_json(directory),
            id="invalid-json",
        ),
        pytest.param(
            lambda directory: directory / "missing.json",
            id="missing-file",
        ),
    ),
)
def test_make_md_summary__fails_for_invalid_input(
    tmp_path: Path, input_factory: Callable[[Path], Path]
) -> None:
    """Invalid input reports should terminate the CLI with an error."""

    input_path = input_factory(tmp_path)
    output_path = tmp_path / "summary.md"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "tools.make_md_summary",
            "--input",
            str(input_path),
            "--output",
            str(output_path),
        ],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 1
    assert result.stdout == ""
    assert output_path.exists() is False
    assert "Input report" in result.stderr


def _create_invalid_json(directory: Path) -> Path:
    path = directory / "report.json"
    path.write_text("{not-json", encoding="utf-8")
    return path
