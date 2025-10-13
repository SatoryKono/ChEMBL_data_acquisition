"""Unit tests for the reporting policy helpers."""

from __future__ import annotations

from typing import Any

import pytest

from tests.helpers.reporting_policy import (
    PIPELINE_SCENARIOS,
    SUCCESS_RATE_THRESHOLD,
    ReportingPolicyState,
)


@pytest.mark.unit
def test_reporting_policy_state__computes_summary_and_success_rate() -> None:
    state = ReportingPolicyState()
    state.register_item("tests/unit/test_demo.py::test_pass", ["csv_loading"])
    state.record_status("tests/unit/test_demo.py::test_pass", "passed")
    state.record_status("tests/unit/test_demo.py::test_skip", "skipped")

    state.finalise(duration_sec=1.234)

    summary = state.build_summary()
    assert summary["total"] == 2
    assert summary["passed"] == 1
    assert summary["skipped"] == 1
    assert summary["success_rate"] == pytest.approx(1.0)
    assert state.missing_scenarios == [
        scenario for scenario in sorted(PIPELINE_SCENARIOS) if scenario != "csv_loading"
    ]


@pytest.mark.unit
def test_reporting_policy_state__flags_missing_scenario_and_threshold_failure() -> None:
    state = ReportingPolicyState(success_threshold=SUCCESS_RATE_THRESHOLD)
    state.register_item("tests/unit/test_demo.py::test_case", ["csv_loading"])
    state.record_status("tests/unit/test_demo.py::test_case", "failed")

    state.finalise(duration_sec=0.5)

    assert "csv_loading" in state.missing_scenarios
    assert state.failure_reasons
    assert any("Success rate" in reason for reason in state.failure_reasons)
    assert any(
        "Missing pipeline scenarios" in reason for reason in state.failure_reasons
    )


@pytest.mark.unit
def test_reporting_policy_state__updates_json_report_payload() -> None:
    state = ReportingPolicyState()
    nodeid = "tests/unit/test_demo.py::test_case"
    state.register_item(nodeid, ["csv_loading"])
    state.record_status(nodeid, "passed")
    state.finalise(duration_sec=2.5)

    report: dict[str, Any] = {"tests": [{"nodeid": nodeid, "status": "passed"}]}
    state.update_json_report(
        report,
        repo="example/repo",
        commit="abcdef1",
        branch="feature/policy",
        python_version="3.11.6",
        pytest_version="7.4.0",
    )

    summary = report["summary"]
    assert summary["total"] == 1
    assert summary["passed"] == 1
    assert summary["success_rate"] == pytest.approx(1.0)

    meta = report["meta"]
    assert meta["repo"] == "example/repo"
    assert meta["commit"] == "abcdef1"
    assert meta["branch"] == "feature/policy"
    assert meta["pipeline_scenarios"]["csv_loading"]["covered"] is True


@pytest.mark.unit
def test_reporting_policy_state__rejects_unknown_scenario() -> None:
    state = ReportingPolicyState()
    with pytest.raises(ValueError):
        state.register_item("tests/test_demo.py::test_case", ["unknown-scenario"])
