from __future__ import annotations

import argparse
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
        gad,
        "write_csv_deterministic",
        lambda df, output, *, key_cols, col_order, chunksize, sort_chunksize, sep, encoding, cfg, **__: output,
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

    def fake_write_csv(
        df,
        output,
        *,
        key_cols,
        col_order,
        chunksize,
        sort_chunksize,
        sep,
        encoding,
        cfg,
        **__,
    ) -> Path:
        captured["col_order"] = list(col_order or [])
        return output

    monkeypatch.setattr(gad, "write_csv_deterministic", fake_write_csv)

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


def test_run_chembl_streams_large_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Ensure large datasets trigger chunked deterministic CSV writes."""

    chunk_size = cfg.io.csv_chunksize
    total_rows = chunk_size * 2 + 5
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n" + "\n".join(str(i) for i in range(total_rows)))

    monkeypatch.setattr(
        io,
        "read_ids",
        lambda *_, **__: (str(i) for i in range(total_rows)),
    )

    def fake_get(ids, cfg, client, chunk_size, timeout, **kwargs):
        data = list(ids)
        return pd.DataFrame(
            {
                "activity_id": data,
                "standard_value": [1.0] * len(data),
                "standard_type": ["IC50"] * len(data),
            }
        )

    monkeypatch.setattr(cl, "get_activities", fake_get)
    monkeypatch.setattr(gad, "normalize_activities", lambda df: df)
    monkeypatch.setattr(gad, "add_pipeline_metadata", lambda df: df)
    monkeypatch.setattr(gad, "compute_activity_bounds", lambda df, cfg: df)
    monkeypatch.setattr(
        gad,
        "apply_activity_annotations",
        lambda df, action_cfg, properties_cfg: df,
    )

    class _Result:
        def __init__(self, data: pd.DataFrame) -> None:
            self.data = data
            self.failure_cases = pd.DataFrame()

    monkeypatch.setattr(
        gad,
        "validate_activities",
        lambda df, return_result: _Result(df),
    )
    monkeypatch.setattr(gad, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gad, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gad, "file_sha256", lambda p: "deadbeef")

    captured: dict[str, object] = {}

    def fake_write_csv(
        df: pd.DataFrame,
        output: Path,
        *,
        key_cols,
        col_order,
        chunksize,
        sort_chunksize,
        sep,
        encoding,
        cfg,
        **__,
    ) -> Path:
        chunk_count = (len(df) + chunksize - 1) // chunksize if chunksize else 1
        captured.update(
            {
                "chunksize": chunksize,
                "sort_chunksize": sort_chunksize,
                "chunk_count": chunk_count,
                "rows": len(df),
            }
        )
        return output

    monkeypatch.setattr(gad, "write_csv_deterministic", fake_write_csv)

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")
    rc = gad.run_chembl(cfg, args)

    assert rc == 0
    assert captured["rows"] == total_rows
    assert captured["chunksize"] == chunk_size
    assert captured["sort_chunksize"] == chunk_size
    assert captured["chunk_count"] > 1
