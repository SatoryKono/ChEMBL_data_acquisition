"""Unit tests for :mod:`library.reporting.test_summary`."""

from __future__ import annotations

import pytest

from library.reporting import test_summary


@pytest.mark.unit
def test_build_summary_markdown__renders_failure_block() -> None:
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

    summary_md = test_summary.build_summary_markdown(report)

    assert "# Test Summary" in summary_md
    assert "- `tests/unit/test_example.py::test_failure` (failed)" in summary_md
    assert summary_md.count("```") == 2
    assert "AssertionError: boom" in summary_md
    assert "line 1" in summary_md and "line 2" in summary_md


@pytest.mark.unit
@pytest.mark.parametrize(
    "report,expected_message",
    [
        ({"summary": {}, "tests": []}, "Summary missing required keys"),
        (
            {
                "summary": {
                    "total": "one",
                    "passed": 0,
                    "failed": 0,
                    "skipped": 0,
                    "xfailed": 0,
                    "xpassed": 0,
                    "error": 0,
                    "success_rate": 1.0,
                },
                "tests": [],
            },
            "must be a non-negative integer",
        ),
        (
            {
                "summary": {
                    "total": 1,
                    "passed": 1,
                    "failed": 0,
                    "skipped": 0,
                    "xfailed": 0,
                    "xpassed": 0,
                    "error": 0,
                    "success_rate": 1.0,
                },
                "tests": {},
            },
            "tests section must be a list",
        ),
    ],
)
def test_build_summary_markdown__validation_errors(
    report: dict[str, object], expected_message: str
) -> None:
    with pytest.raises(ValueError, match=expected_message):
        test_summary.build_summary_markdown(report)
