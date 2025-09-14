"""Tests for validate_and_write utility."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from pytest import MonkeyPatch

from library.config import Config
from scripts import get_target_data as gtd


def test_validate_and_write(
    tmp_path: Path, cfg: Config, monkeypatch: MonkeyPatch
) -> None:
    df = pd.DataFrame({"target_chembl_id": ["C1"]})
    out = tmp_path / "out.csv"
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *a, **k: None)
    exit_code = gtd.validate_and_write(df, out, cfg)
    assert exit_code == 0
    assert out.exists()
