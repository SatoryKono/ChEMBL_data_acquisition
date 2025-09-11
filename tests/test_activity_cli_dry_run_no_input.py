from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

from library.cli import LoggerConfig, configure_logger
from library.config import Config
from scripts import get_activity_data as gad


def test_dry_run_no_input(tmp_path: Path) -> None:
    buf = io.StringIO()
    configure_logger(LoggerConfig(stream=buf))
    args = argparse.Namespace(
        input_csv=Path("missing.csv"),
        output_csv=None,
    )
    cfg = Config()
    cfg.activity.limit = 5
    cfg.activity.dry_run = True
    rc = gad.run_chembl(cfg, args)
    assert rc == 0
    lines = buf.getvalue().splitlines()
    assert lines
    record = json.loads(lines[-1])
    assert record.get("event") == "dry_run"
    assert "dry run selected" in record.get("msg", "")
    assert "5 identifiers" in record.get("msg", "")
