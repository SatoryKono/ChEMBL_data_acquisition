"""Tests for :mod:`library.postprocessing.common.logging`."""

from __future__ import annotations

from library.postprocessing.common.logging import (
    PipelineRunMetrics,
    build_report_payload,
)


def _make_metrics() -> PipelineRunMetrics:
    metrics = PipelineRunMetrics(
        pipeline_version="1.0",
        started_at="2024-01-01T00:00:00+00:00",
        input_rows=10,
        input_columns=5,
    )
    metrics.output_rows = 8
    metrics.output_columns = 4
    metrics.completed_at = "2024-01-01T00:01:00+00:00"
    metrics.duration_s = 60.0
    return metrics


def test_build_report_payload__includes_legacy_output_key():
    metrics = _make_metrics()

    payload = build_report_payload(
        table="activities",
        metrics=metrics,
        output_path="legacy.csv",
    )

    assert payload["output_path"] == "legacy.csv"
    assert payload["output_postprocessed"] == "legacy.csv"


def test_build_report_payload__prefers_explicit_postprocessed_value():
    metrics = _make_metrics()

    payload = build_report_payload(
        table="activities",
        metrics=metrics,
        output_path="legacy.csv",
        output_postprocessed="postprocessed.csv",
    )

    assert payload["output_path"] == "legacy.csv"
    assert payload["output_postprocessed"] == "postprocessed.csv"


def test_build_report_payload__omits_legacy_key_when_not_available():
    metrics = _make_metrics()

    payload = build_report_payload(
        table="activities",
        metrics=metrics,
        output_path=None,
        output_postprocessed="postprocessed.csv",
    )

    assert "output_path" not in payload
    assert payload["output_postprocessed"] == "postprocessed.csv"
