"""Tests for helper utilities in :mod:`scripts.get_target_data`."""

from pathlib import Path

import pandas as pd

from library.config import Config
from scripts.get_target_data import _first_token, _pipe_merge, _save_snapshot


def test_pipe_merge_deduplicates_and_sorts() -> None:
    values = ["b|a", "C|a", None, "", "b"]
    assert _pipe_merge(values) == "C|a|b"


def test_pipe_merge_handles_empty() -> None:
    assert _pipe_merge([None, "", " "]) == ""


def test_first_token_extracts_head() -> None:
    assert _first_token("a|b|c") == "a"
    assert _first_token(None) == ""
    assert _first_token("") == ""


def test_save_snapshot_creates_unique_files(tmp_path: Path, cfg: Config) -> None:
    df = pd.DataFrame({"a": [1]})
    base = tmp_path / "out.csv"
    _save_snapshot(df, base, "step", cfg)
    _save_snapshot(df, base, "step", cfg)
    files = sorted(tmp_path.glob("out_step_*.csv"))
    assert len(files) == 2
