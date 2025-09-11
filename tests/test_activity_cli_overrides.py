from __future__ import annotations

from pathlib import Path

import get_activity_data as gad
from library import chembl_library as cl, io
import pandas as pd


def _create_config(tmp_path: Path) -> Path:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "jobs:\n  chunk_size: 10\n"
        "io:\n  csv_sep: '|'\n  csv_encoding: iso-8859-1\n"
        "log:\n  level: INFO\n"
        "api:\n  timeout_read: 30\n"
        "resources:\n"
        "  dictionary_dir: dictionary\n"
        "  iuphar_target_csv: dictionary/_IUPHAR/_IUPHAR_target.csv\n"
        "  iuphar_family_csv: dictionary/_IUPHAR/_IUPHAR_family.csv\n"
        "  uniprot_data_dir: uniprot\n"
    )
    return cfg


def _run(
    tmp_path: Path, monkeypatch, extra: list[str]
) -> tuple[int, dict[str, object]]:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n")
    config_path = _create_config(tmp_path)
    called: dict[str, object] = {}
    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["1"])

    def fake_get(ids, cfg, chunk_size, timeout):
        called["chunk_size"] = chunk_size
        return pd.DataFrame({"activity_id": ids})

    def fake_write(df, output, cfg, sep, encoding):
        called["sep"] = sep
        called["encoding"] = encoding

    monkeypatch.setattr(cl, "get_activities", fake_get)
    monkeypatch.setattr(io, "write_csv", fake_write)
    monkeypatch.setattr(gad, "analyze_table_quality", lambda df, table_name: None)
    rc = gad.main(["--config", str(config_path), "--input", str(input_csv), *extra])
    return rc, called


def test_default_config_used(tmp_path: Path, monkeypatch) -> None:
    rc, called = _run(tmp_path, monkeypatch, [])
    assert rc == 0
    assert called["chunk_size"] == 10
    assert called["sep"] == "|"
    assert called["encoding"] == "iso-8859-1"


def test_cli_overrides(tmp_path: Path, monkeypatch) -> None:
    rc, called = _run(
        tmp_path,
        monkeypatch,
        ["--chunk-size", "3", "--sep", ";", "--encoding", "latin1"],
    )
    assert rc == 0
    assert called["chunk_size"] == 3
    assert called["sep"] == ";"
    assert called["encoding"] == "latin1"
