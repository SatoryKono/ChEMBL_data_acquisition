"""Unit tests for helpers in :mod:`library.cli.entrypoints.activity`."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library.cli.entrypoints import activity


@pytest.mark.unit
@pytest.mark.parametrize(
    ("dtype", "expected_dtype", "values"),
    [
        ("Float64", "Float64", [1.5, None, 2.0]),
        ("Int64", "Int64", [1, None, 3]),
        ("boolean", "boolean", [True, None, False]),
        ("string", "string", ["foo", None, "bar"]),
    ],
)
def test_coerce_series_dtype__extension_roundtrip(
    dtype: str, expected_dtype: str, values: list[object]
) -> None:
    series = pd.Series(values)

    result = activity._coerce_series_dtype(series, dtype)

    assert str(result.dtype) == expected_dtype
    expected_series = pd.Series(values, dtype=expected_dtype)
    pd.testing.assert_series_equal(result, expected_series)


@pytest.mark.unit
@pytest.mark.parametrize(
    "dtype",
    ["category", "invalid-type"],
)
def test_coerce_series_dtype__fallback_to_string(dtype: str) -> None:
    series = pd.Series(["1", "2", None])

    result = activity._coerce_series_dtype(series, dtype)

    assert result.dtype == object
    assert result.tolist() == ["1", "2", "None"]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("value", "default", "expected"),
    [
        (5, 0, 5),
        ("7", 0, 7),
        (None, 3, 3),
        ("not-a-number", 2, 2),
        (["1"], 4, 4),
    ],
)
def test_safe_int__conversion_cases(value: object, default: int, expected: int) -> None:
    assert activity._safe_int(value, default) == expected


class _StubLogger:
    def __init__(self) -> None:
        self.calls: list[tuple[str, dict[str, Any]]] = []
        self.events: list[tuple[str, str, dict[str, Any]]] = []

    def info(self, event: str, **payload: Any) -> None:
        self.calls.append((event, payload))


@pytest.mark.unit
def test_emit_completion_message__basic_payload(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    logger = _StubLogger()
    monkeypatch.setattr(activity, "logger", logger)

    activity._emit_completion_message(
        output_path=Path("output.csv"),
        processed_rows=15,
        duration_s=1.25,
        mode="fetch",
    )

    assert len(logger.calls) == 1
    event, payload = logger.calls[0]
    assert event == "activity_pipeline_completion"
    assert payload == {
        "output": "output.csv",
        "rows": 15,
        "duration_s": 1.25,
        "mode": "fetch",
        "processed_rows": 15,
    }


@pytest.mark.unit
def test_emit_completion_message__skip_existing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    logger = _StubLogger()
    monkeypatch.setattr(activity, "logger", logger)

    activity._emit_completion_message(
        output_path=Path("existing.csv"),
        processed_rows=None,
        duration_s=0.5,
        mode="skip_existing",
    )

    assert len(logger.calls) == 2
    skip_event, skip_payload = logger.calls[0]
    assert skip_event == "pipeline_skip_existing"
    assert skip_payload == {"output": "existing.csv"}

    completion_event, completion_payload = logger.calls[1]
    assert completion_event == "activity_pipeline_completion"
    assert completion_payload == {
        "output": "existing.csv",
        "rows": 0,
        "duration_s": 0.5,
        "mode": "skip_existing",
    }

    assert logger.events == [
        ("info", "pipeline_skip_existing", {"output": "existing.csv"}),
        (
            "info",
            "activity_pipeline_completion",
            {
                "output": "existing.csv",
                "rows": 0,
                "duration_s": 0.5,
                "mode": "skip_existing",
            },
        ),
    ]


@pytest.mark.unit
def test_emit_completion_message__streamed_metrics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    logger = _StubLogger()
    monkeypatch.setattr(activity, "logger", logger)

    class ReprObject:
        def __str__(self) -> str:  # pragma: no cover - trivial
            return "repr-object"

    activity._emit_completion_message(
        output_path=None,
        processed_rows=None,
        duration_s=2.0,
        mode="stream",
        streamed_metrics={
            "rows": 9.7,
            "null_fraction": 0.125,
            "nan_value": float("nan"),
            "integral_value": 4,
            "text": "value",
            "flag": True,
            "other": ReprObject(),
        },
    )

    assert len(logger.calls) == 1
    event, payload = logger.calls[0]
    assert event == "activity_pipeline_completion"
    assert payload["output"] is None
    assert payload["rows"] == 9
    assert payload["duration_s"] == 2.0
    assert payload["mode"] == "stream"
    assert payload["null_fraction"] == pytest.approx(0.125)

    metrics = payload["streamed_metrics"]
    assert metrics["rows"] == pytest.approx(9.7)
    assert metrics["null_fraction"] == pytest.approx(0.125)
    assert metrics["nan_value"] is None
    assert metrics["integral_value"] == 4
    assert metrics["text"] == "value"
    assert metrics["flag"] == 1
    assert metrics["other"] == "repr-object"
