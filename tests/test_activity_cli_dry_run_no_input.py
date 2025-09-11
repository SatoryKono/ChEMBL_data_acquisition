from __future__ import annotations

from pathlib import Path
import io
import json
import argparse

from library.cli import LoggerConfig, configure_logger
import get_activity_data as gad
from library.config import Config

 
def test_dry_run_no_input(tmp_path: Path) -> None:
    buf = io.StringIO()
    configure_logger(LoggerConfig(stream=buf))
    args = argparse.Namespace(
        input_csv=Path("missing.csv"),
        output_csv=None,
        column="activity_id",
        sep=",",
        encoding="utf8",
        chunk_size=5,
        timeout=30.0,
        limit=5,
        dry_run=True,
 
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
    rc = gad.run_chembl(Config(), args)
    assert rc == 0
    lines = buf.getvalue().splitlines()
    assert lines
    record = json.loads(lines[-1])
    assert "dry run selected" in record.get("msg", "")
    assert "5 identifiers" in record.get("msg", "")
