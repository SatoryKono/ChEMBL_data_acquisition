"""Unit tests for helpers in :mod:`library.cli.entrypoints.activity`."""

from __future__ import annotations

import argparse
import csv
import importlib
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library.cli.entrypoints import activity
from library.config import Config


@pytest.mark.unit
def test_activity_runner_registration__commands_module(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    commands_module = importlib.import_module("library.cli.commands.get_activity_data")
    runner_module = importlib.import_module("library.pipelines.activity.runner")

    def _stub_runner(cfg, args):  # pragma: no cover - exercised via resolve
        return 0

    def _stub_emit(**_kwargs):  # pragma: no cover - exercised via resolve
        return None

    monkeypatch.setattr(commands_module, "run_chembl", _stub_runner)
    monkeypatch.setattr(commands_module, "_emit_completion_message", _stub_emit)

    reloaded_runner = importlib.reload(runner_module)

    try:
        runner, emit_completion = reloaded_runner.resolve_activity_pipeline_hooks()
        assert runner is _stub_runner
        assert emit_completion is _stub_emit
    finally:
        monkeypatch.undo()
        importlib.reload(reloaded_runner)


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
        self.events.append(("info", event, dict(payload)))


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
        "output_postprocessed": "output.csv",
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

    assert logger.calls == [
        ("pipeline_skip_existing", {"output_postprocessed": "existing.csv"}),
        (
            "activity_pipeline_completion",
            {
                "output": "existing.csv",
                "rows": 0,
                "duration_s": 0.5,
                "mode": "skip_existing",
            },
        ),
    ]
    assert (
        "info",
        "pipeline_skip_existing",
        {"output_postprocessed": "existing.csv"},
    ) in logger.events
    assert any(
        event == "activity_pipeline_completion"
        and payload.get("mode") == "skip_existing"
        for _, event, payload in logger.events
    )


def test_run_activity_pipeline__propagates_emit_legacy(
    cfg: Config, tmp_path: Path
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    options = activity.ActivityCommandOptions(
        input_csv=input_csv,
        output_csv=output_csv,
        final_output=output_csv,
        emit_legacy_artifacts=True,
    )

    recorded: dict[str, Any] = {}

    def fake_runner(_cfg: Config, args: argparse.Namespace) -> int:
        recorded["emit"] = getattr(args, "emit_legacy_artifacts", False)
        return 0

    exit_code = activity.run_activity_pipeline(
        cfg,
        options,
        runner=fake_runner,
        emit_completion_message=lambda **_: None,
    )

    assert exit_code == 0
    assert recorded["emit"] is True


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
    assert payload["output_postprocessed"] is None
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


@pytest.mark.unit
def test_resolve_completion_rows__prefers_pipeline_stats_when_summary_zero() -> None:
    result = activity._resolve_completion_rows(
        processed_ids=10,
        summary_snapshot={"rows": 0},
        pipeline_stats={"rows_kept": 10, "rows_total": 10},
    )

    assert result == 10


@pytest.mark.unit
def test_resolve_completion_rows__falls_back_to_summary_when_stats_missing() -> None:
    result = activity._resolve_completion_rows(
        processed_ids=8,
        summary_snapshot={"rows": 5},
        pipeline_stats=None,
    )

    assert result == 5


@pytest.mark.unit
def test_resolve_completion_rows__falls_back_to_processed_when_metrics_absent() -> None:
    result = activity._resolve_completion_rows(
        processed_ids=7,
        summary_snapshot=None,
        pipeline_stats=None,
    )

    assert result == 7


@pytest.mark.unit
def test_load_assay_src_lookup__coerces_numeric_values(tmp_path: Path) -> None:
    dictionary_dir = tmp_path / "dictionary"
    assay_dir = dictionary_dir / "_assay"
    assay_dir.mkdir(parents=True)
    csv_path = assay_dir / "assay.csv"

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["assay_chembl_id", "src_assay_id"])
        writer.writerow(["CHEMBL_INT", 357280])
        writer.writerow(["CHEMBL_FLOAT", 357280.0])
        writer.writerow(["CHEMBL_STR", "  001  "])
        writer.writerow(["CHEMBL_EMPTY", "   "])
        writer.writerow(["", "12345"])

    lookup = activity._load_assay_src_lookup(dictionary_dir)

    assert lookup["CHEMBL_INT"] == "357280"
    assert lookup["CHEMBL_FLOAT"] == "357280"
    assert lookup["CHEMBL_STR"] == "001"
    assert "CHEMBL_EMPTY" not in lookup
    assert "" not in lookup


@pytest.mark.unit
def test_derive_standard_output_labels__handles_hidden_tmp_suffix() -> None:
    table, date = activity._derive_standard_output_labels(
        Path(".output.activities_20240101.csv.tmp")
    )

    assert table == "activity"
    assert date == "20240101"


@pytest.mark.unit
def test_derive_standard_output_labels__strips_output_prefix() -> None:
    table, date = activity._derive_standard_output_labels(
        Path("output.documents_20231231.csv")
    )

    assert table == "documents"
    assert date == "20231231"


@pytest.mark.unit
def test_derive_standard_output_labels__deduplicates_output_prefix_chain() -> None:
    table, date = activity._derive_standard_output_labels(
        Path("output..output.activities_20251011.csv.tmp")
    )

    assert table == "activity"
    assert date == "20251011"


@pytest.mark.unit
def test_derive_standard_output_labels__fallbacks_to_current_date(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(activity, "_current_date_token", lambda: "19990101")

    table, date = activity._derive_standard_output_labels(Path("custom_export.csv"))

    assert table == "custom"
    assert date == "19990101"
