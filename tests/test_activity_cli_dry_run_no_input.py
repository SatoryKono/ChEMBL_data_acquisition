from __future__ import annotations

from pathlib import Path

import get_activity_data as gad


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


def test_dry_run_no_input(tmp_path: Path, capfd) -> None:
    config_path = _create_config(tmp_path)
    rc = gad.main(["--config", str(config_path), "--dry-run", "--limit", "5"])
    assert rc == 0
    stderr = capfd.readouterr().err
    assert "dry run selected" in stderr
    assert "5 identifiers" in stderr
