from __future__ import annotations

from pathlib import Path

import yaml

from scripts.table_quality_main import main


def test_print_config_cli(tmp_path: Path, monkeypatch, capsys) -> None:
    csv_path = tmp_path / "data.csv"
    csv_path.write_text("chembl_id\nCHEMBL1\n")
    config = Path("tests/fixtures/config.min.yaml")
    monkeypatch.setenv("CHEMBL_DA__API__RPS", "7")
    monkeypatch.setenv("CHEMBL_DA__LOG__LEVEL", "INFO")
    rc = main(
        [
            "--config",
            str(config),
            "--input",
            str(csv_path),
            "--table-name",
            "demo",
            "--log-level",
            "DEBUG",
            "--output",
            str(tmp_path / "out"),
            "--print-config",
        ]
    )
    assert rc == 0
    out = capsys.readouterr().out
    data = yaml.safe_load(out)
    assert data["api"]["rps"] == 7
    assert data["log"]["level"] == "DEBUG"
    assert not (tmp_path / "out").exists()
