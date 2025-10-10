"""Unit tests for :mod:`library.reporting.run_manifest`."""

from __future__ import annotations

from pathlib import Path

import pytest

from library.reporting import run_manifest


@pytest.mark.unit
def test_finalise_csv_output__returns_report_with_expected_stats(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "result.csv"
    csv_path.write_text("id\n1\n2\n", encoding="utf-8")

    report = run_manifest.finalise_csv_output(
        csv_path=csv_path,
        rows_total=2,
        rows_kept=2,
        command="unit-test",
        config={"system": {"flag": True}},
        inputs={"input_csv": "source.csv"},
        schema="TestSchema",
        stats_extra={"custom_metric": 3},
        quality_summary={"rows_total": 2},
    )

    stats = report.stats
    assert stats["rows_total"] == 2
    assert stats["rows_kept"] == 2
    assert stats["rows_dropped"] == 0
    assert stats["custom_metric"] == 3

    assert report.meta_path.name == "result.csv.meta.yaml"
    assert report.meta_sha256

    assert report.quality_path is not None
    assert report.quality_path.name == "result.csv.quality.json"
    assert report.quality_sha256

    manifest_payload = report.as_manifest_payload()
    assert manifest_payload["stats"]["rows_total"] == 2
    assert manifest_payload["stats"]["rows_kept"] == 2
    assert manifest_payload["stats"]["rows_dropped"] == 0
    assert manifest_payload["stats"]["custom_metric"] == 3
    assert manifest_payload["output_path"].endswith("result.csv")
    assert manifest_payload["meta"]["path"].endswith("result.csv.meta.yaml")

    entry = {"output": {"exists": True}}
    run_manifest.merge_run_output(entry, report)
    assert entry["output"]["exists"] is True
    assert entry["output"]["meta_path"].endswith("result.csv.meta.yaml")
    assert entry["output"]["quality_path"].endswith("result.csv.quality.json")
    assert entry["stats"]["rows_total"] == 2
    assert entry["stats"]["rows_kept"] == 2
    assert entry["stats"]["rows_dropped"] == 0

    loaded = run_manifest.load_output_report(csv_path)
    assert loaded is not None
    assert loaded.stats["rows_total"] == 2
    assert loaded.meta_path == report.meta_path


@pytest.mark.unit
def test_finalise_csv_output__quality_builder_invoked(tmp_path: Path) -> None:
    csv_path = tmp_path / "data.csv"
    csv_path.write_text("id\n1\n", encoding="utf-8")

    calls: list[dict[str, int]] = []

    def _builder(summary: dict[str, int]) -> dict[str, int]:
        calls.append(summary)
        return {"value": summary["value"]}

    report = run_manifest.finalise_csv_output(
        csv_path=csv_path,
        rows_total=1,
        rows_kept=1,
        command="builder",
        config={},
        inputs={},
        schema="BuilderSchema",
        quality_summary={"value": 5},
        quality_builder=_builder,
    )

    assert calls == [{"value": 5}]
    assert report.quality_path is not None
    assert report.quality_path.name == "data.csv.quality.json"


@pytest.mark.unit
def test_finalise_csv_output__quality_error(tmp_path: Path) -> None:
    csv_path = tmp_path / "invalid.csv"
    csv_path.write_text("id\n1\n", encoding="utf-8")

    with pytest.raises(run_manifest.QualityReportError) as excinfo:
        run_manifest.finalise_csv_output(
            csv_path=csv_path,
            rows_total=1,
            rows_kept=1,
            command="quality",
            config={},
            inputs={},
            schema="ErrorSchema",
            quality_summary={"bad": {1, 2}},
        )

    assert excinfo.value.path == csv_path.with_suffix(".quality.json")


@pytest.mark.unit
def test_finalise_csv_output__quality_analysis_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    csv_path = tmp_path / "analysis.csv"
    csv_path.write_text("id\n1\n", encoding="utf-8")

    def _raise(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("quality boom")

    monkeypatch.setattr(run_manifest, "analyze_table_quality", _raise)

    with pytest.raises(run_manifest.QualityAnalysisError):
        run_manifest.finalise_csv_output(
            csv_path=csv_path,
            rows_total=1,
            rows_kept=1,
            command="analysis",
            config={},
            inputs={},
            schema="AnalysisSchema",
            quality_summary={},
            quality_profiler=object(),
            quality_config={"enable": True},
        )


@pytest.mark.unit
def test_load_output_report__missing(tmp_path: Path) -> None:
    csv_path = tmp_path / "missing.csv"
    assert run_manifest.load_output_report(csv_path) is None
