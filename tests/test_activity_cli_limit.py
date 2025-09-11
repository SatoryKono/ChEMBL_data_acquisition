from __future__ import annotations

from pathlib import Path

import pandas as pd

import get_activity_data as gad
from library import chembl_library as cl, io


def _create_config(tmp_path: Path) -> Path:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "jobs:\n  chunk_size: 10\n"
        "io:\n  csv_sep: ','\n  csv_encoding: utf8\n"
        "log:\n  level: INFO\n"
        "api:\n  timeout_read: 30\n"
        "resources:\n"
        "  dictionary_dir: dictionary\n"
        "  iuphar_target_csv: dictionary/_IUPHAR/_IUPHAR_target.csv\n"
        "  iuphar_family_csv: dictionary/_IUPHAR/_IUPHAR_family.csv\n"
        "  uniprot_data_dir: uniprot\n"
        "  organism_csv: dictionary/organism.csv\n"
        "  status_csv: dictionary/status.csv\n"
        "  targets_type_csv: dictionary/targets_type.csv\n"
    )
    return cfg


def test_activity_cli_limit(tmp_path: Path, monkeypatch) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n2\n3\n")
    config_path = _create_config(tmp_path)

    called: dict[str, object] = {"ids": None, "written": False}
    monkeypatch.setattr(io, "read_ids", lambda *a, **k: iter(["1", "2", "3"]))

    def fake_get(ids, cfg, chunk_size, timeout):
        data = list(ids)
        called["ids"] = data
        return pd.DataFrame({"activity_id": data})

    def fake_write(df, output, cfg, sep, encoding):
        called["written"] = True

    monkeypatch.setattr(cl, "get_activities", fake_get)
    monkeypatch.setattr(io, "write_csv", fake_write)
    monkeypatch.setattr(gad, "analyze_table_quality", lambda df, table_name: None)

    rc = gad.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--limit",
            "2",
        ]
    )
    assert rc == 0
    assert called["ids"] == ["1", "2"]

    called["ids"] = None
    called["written"] = False
    rc = gad.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--limit",
            "2",
            "--dry-run",
        ]
    )
    assert rc == 0
    assert called["ids"] is None
    assert called["written"] is False
