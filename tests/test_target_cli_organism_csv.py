from __future__ import annotations

from pathlib import Path

import get_target_data as gtd


def test_organism_csv_default_from_config(tmp_path, monkeypatch):
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        "jobs:\n  chunk_size: 10\n"
        "io:\n  csv_sep: ','\n  csv_encoding: utf8\n"
        "log:\n  level: INFO\n"
        "api:\n  timeout_read: 30\n"
        "resources:\n"
        "  dictionary_dir: dictionary\n"
        "  iuphar_target_csv: dictionary/_IUPHAR/_IUPHAR_target.csv\n"
        "  iuphar_family_csv: dictionary/_IUPHAR/_IUPHAR_family.csv\n"
        "  organism_csv: dictionary/organism.csv\n"
    )
    captured: dict[str, Path] = {}

    def fake_run_all(cfg, args):  # type: ignore[unused-argument]
        captured["organism_csv"] = args.organism_csv
        return 0

    monkeypatch.setattr(gtd, "run_all", fake_run_all)
    rc = gtd.main(["all", "--config", str(cfg_path)])
    assert rc == 0
    assert captured["organism_csv"] == Path("dictionary/organism.csv")
