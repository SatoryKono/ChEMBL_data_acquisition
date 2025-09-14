"""Tests for missing UniProt column handling in :mod:`scripts.get_target_data`."""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

import pandas as pd
import pytest

from library.cli import LoggerConfig, configure_logger
from library.config import Config
from scripts import get_target_data as gtd


def test_run_all_missing_uniprot_column(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Run fails when the expected UniProt column is absent."""
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n")

    def fake_run_chembl(cfg: Config, args: argparse.Namespace) -> int:
        pd.DataFrame({"target_chembl_id": ["CHEMBL1"]}).to_csv(
            args.output_csv, index=False
        )
        return 0

    monkeypatch.setattr(gtd, "run_chembl", fake_run_chembl)
    monkeypatch.setattr(
        gtd, "run_uniprot", lambda cfg, args: pytest.fail("unexpected call")
    )
    monkeypatch.setattr(
        gtd, "run_iuphar", lambda cfg, args: pytest.fail("unexpected call")
    )

    buffer = io.StringIO()
    configure_logger(LoggerConfig(level="ERROR", stream=buffer))

    cfg = Config()
    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    rc = gtd.run_all(cfg, args)
    assert rc == 1
    log_output = buffer.getvalue()
    configure_logger(LoggerConfig(stream=sys.stdout))
    assert "missing column" in log_output and "uniprot_id" in log_output
