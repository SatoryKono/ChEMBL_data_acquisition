from __future__ import annotations

import argparse
import types
from pathlib import Path

import pandas as pd
import pytest

from library import chembl_library as cl
from library import io
from library.config import Config
from schemas import ActivitiesSchema
from scripts import get_activity_data as gad


def test_run_chembl_respects_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n2\n3\n")

    cfg.activity.limit = 2
    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    read_ids: list[str] = []

    def fake_read_ids(*_, **__):
        for item in ["1", "2", "3"]:
            read_ids.append(item)
            yield item

    monkeypatch.setattr(io, "read_ids", fake_read_ids)

    captured: dict[str, object] = {}

    def fake_get(ids, cfg, client, chunk_size, timeout, **kwargs):
        data = list(ids)
        captured["ids"] = data
        captured["extra_columns"] = kwargs.get("extra_columns")
        return pd.DataFrame({"activity_id": data})

    monkeypatch.setattr(cl, "get_activities", fake_get)

    monkeypatch.setattr(
        io,
        "write_csv",
        lambda df, output, *, cfg, sep=None, encoding=None, **__: output,
    )
    monkeypatch.setattr(gad, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gad, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gad, "file_sha256", lambda p: "deadbeef")

    rc = gad.run_chembl(cfg, args)
    assert rc == 0
    assert captured["ids"] == ["1", "2"]
    assert read_ids == ["1", "2"]
    assert captured["extra_columns"] == ["action_type"]


def test_run_chembl_limit_dry_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n2\n3\n")

    cfg.activity.limit = 2
    cfg.activity.dry_run = True
    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    called = {"read": False}

    def fake_read_ids(*_, **__):
        called["read"] = True
        return iter([])

    monkeypatch.setattr(io, "read_ids", fake_read_ids)

    rc = gad.run_chembl(cfg, args)
    assert rc == 0
    assert called["read"] is False


def test_run_chembl_column_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Ensure schema columns precede alphabetically sorted extras."""
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["1"]))

    df = pd.DataFrame(
        [
            {
                "standard_value": 1.0,
                "assay_description": "desc",
                "activity_id": "1",
                "molecule_chembl_id": "CHEMBL1",
                "bao_format": "A",
                "standard_type": "IC50",
                "target_id": "T1",
                "assay_chembl_id": "A1",
            }
        ]
    )

    monkeypatch.setattr(cl, "get_activities", lambda *_, **__: df)
    monkeypatch.setattr(gad, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gad, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gad, "file_sha256", lambda p: "deadbeef")

    captured: dict[str, list[str]] = {}

    def fake_write_csv(df, output, *, cfg, key_cols=None, col_order=None, **__) -> Path:
        captured["col_order"] = list(col_order or [])
        return output

    monkeypatch.setattr(io, "write_csv", fake_write_csv)

    rc = gad.run_chembl(cfg, args)
    assert rc == 0

    schema_cols = list(ActivitiesSchema.columns)
    available = set(df.columns) | {
        "pipeline_version",
        "timestamp_utc",
        "lower_value",
        "upper_value",
        "activity_properties",
        "action_type",
        "properties_hash",
    }
    expected_head = [c for c in schema_cols if c in available]
    expected_tail = sorted(available - set(schema_cols))
    assert captured["col_order"] == expected_head + expected_tail


def test_run_chembl_skip_quality_omits_reports(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n")
    output_csv = tmp_path / "out.csv"
    args = argparse.Namespace(
        input_csv=input_csv, output_csv=output_csv, enable_quality=False
    )

    cfg.system.reports.enable_quality = True

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["1"]))

    df = pd.DataFrame(
        [
            {
                "activity_id": "1",
                "assay_chembl_id": "A1",
                "molecule_chembl_id": "CHEMBL1",
                "standard_type": "IC50",
                "standard_value": 10.0,
            }
        ]
    )
    monkeypatch.setattr(cl, "get_activities", lambda *_, **__: df)
    monkeypatch.setattr(gad, "apply_activity_annotations", lambda frame, **__: frame)
    monkeypatch.setattr(gad, "compute_activity_bounds", lambda frame, _: frame)
    monkeypatch.setattr(gad, "add_pipeline_metadata", lambda frame: frame)
    def fake_write_csv(frame: pd.DataFrame, output: Path | str, *, cfg, **__) -> Path:
        path = Path(output)
        frame.to_csv(path, index=False)
        return path

    monkeypatch.setattr(io, "write_csv", fake_write_csv)
    monkeypatch.setattr(gad, "write_meta_yaml", lambda **__: None)
    monkeypatch.setattr(gad, "file_sha256", lambda _: "deadbeef")
    monkeypatch.setattr(
        gad,
        "validate_activities",
        lambda frame, return_result: types.SimpleNamespace(
            data=frame, failure_cases=pd.DataFrame()
        ),
    )

    rc = gad.run_chembl(cfg, args)
    assert rc == 0

    base = output_csv.with_suffix("")
    expected_reports = [
        base.with_name(f"{base.name}_quality_report_table.csv"),
        base.with_name(f"{base.name}_data_correlation_report_table.csv"),
    ]
    for report in expected_reports:
        assert not report.exists()
