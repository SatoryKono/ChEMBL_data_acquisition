from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pytest

from library import chembl_library as cl
from library import io
from library.config import Config
from scripts import get_activity_data as gad


def test_run_chembl_respects_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n2\n3\n")

    cfg = Config()
    cfg.activity.limit = 2
    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    read_ids: list[str] = []

    def fake_read_ids(*_, **__):
        for item in ["1", "2", "3"]:
            read_ids.append(item)
            yield item

    monkeypatch.setattr(io, "read_ids", fake_read_ids)

    captured: dict[str, list[str]] = {}

    def fake_get(ids, cfg, client, chunk_size, timeout):
        data = list(ids)
        captured["ids"] = data
        return pd.DataFrame({"activity_id": data})

    monkeypatch.setattr(cl, "get_activities", fake_get)

    def fake_write_csv_det(df, path, col_order, key_cols):
        path.write_text("id\n1\n")
        return path

    monkeypatch.setattr(gad, "write_csv_deterministic", fake_write_csv_det)
    monkeypatch.setattr(io, "write_csv", lambda df, output, cfg, sep, encoding: None)
    monkeypatch.setattr(gad, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gad, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gad, "file_sha256", lambda p: "deadbeef")

    rc = gad.run_chembl(cfg, args)
    assert rc == 0
    assert captured["ids"] == ["1", "2"]
    assert read_ids == ["1", "2"]


def test_run_chembl_limit_dry_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n2\n3\n")

    cfg = Config()
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
