from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pytest

from library import chembl_library as cl
from library import io
from library.config import Config
from library.processing import activity as activity_processing
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
    monkeypatch.setattr(activity_processing, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(activity_processing, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(activity_processing, "file_sha256", lambda p: "deadbeef")

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
    monkeypatch.setattr(activity_processing, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(activity_processing, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(activity_processing, "file_sha256", lambda p: "deadbeef")

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


def test_main_delegates_to_processing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n")
    output_csv = tmp_path / "out.csv"

    cfg_copy = cfg.model_copy(deep=True)
    cfg_copy.activity.dry_run = False
    cfg_copy.activity.limit = None

    monkeypatch.setattr(gad.cli, "apply_config_overrides", lambda *_, **__: cfg_copy)
    monkeypatch.setattr(gad, "ensure_dirs", lambda _cfg: None)

    called: dict[str, object] = {}

    def fake_fetch(cfg_arg, path, *, limit):
        called["fetch"] = {
            "cfg": cfg_arg,
            "path": path,
            "limit": limit,
        }
        frame_df = pd.DataFrame({"activity_id": ["1"]})
        return activity_processing.ActivityFrameResult(
            data=frame_df,
            column_order=["activity_id"],
            rows_total=1,
        )

    monkeypatch.setattr(activity_processing, "fetch_activity_frame", fake_fetch)
    monkeypatch.setattr(gad, "fetch_activity_frame", fake_fetch)

    def fake_write(frame, output, cfg_arg, *, command, input_csv):
        called["write"] = {
            "frame": frame,
            "output": output,
            "cfg": cfg_arg,
            "command": command,
            "input_csv": input_csv,
        }
        return activity_processing.ActivityWriteResult(
            csv_path=output,
            rows_kept=1,
            rows_dropped=0,
            exit_code=0,
        )

    monkeypatch.setattr(activity_processing, "write_activity_outputs", fake_write)
    monkeypatch.setattr(gad, "write_activity_outputs", fake_write)

    exit_code = gad.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--config",
            "ignored.yaml",
        ]
    )

    assert exit_code == 0
    assert called["fetch"]["path"] == input_csv
    assert called["fetch"]["cfg"] is cfg_copy
    assert called["write"]["output"] == output_csv
    assert called["write"]["cfg"] is cfg_copy
