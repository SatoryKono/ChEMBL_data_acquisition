from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import get_activity_data as gad
from library.cli import LoggerConfig, configure_logger
from library.config import Config


def test_dry_run_no_input(tmp_path: Path) -> None:
    """Dry run should not fail when the input file is missing."""
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
    )
    rc = gad.run_chembl(Config(), args)
    assert rc == 0
    lines = buf.getvalue().splitlines()
    assert lines
    record = json.loads(lines[-1])
    assert "dry run selected" in record.get("msg", "")
    assert "5 identifiers" in record.get("msg", "")
