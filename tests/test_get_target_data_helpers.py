"""Tests for helper utilities in :mod:`scripts.get_target_data`."""

from pathlib import Path

import pandas as pd
import pytest

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
    first = _save_snapshot(df, base, "step", cfg)
    second = _save_snapshot(df, base, "step", cfg)
    files = sorted(tmp_path.glob("out_step_*.csv"))
    assert len(files) == 2
    for path in (first, second):
        meta = path.with_name(path.name + ".meta.yaml")
        assert meta.exists()


def test_save_snapshot_respects_io_settings(tmp_path: Path, cfg: Config) -> None:
    df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "pref_name": ["café"],
        }
    )
    base = tmp_path / "snapshot.csv"
    cfg.io.csv_sep = ";"
    cfg.io.csv_encoding = "latin1"
    path = _save_snapshot(df, base, "stage", cfg)

    with pytest.raises(UnicodeDecodeError):
        path.read_text(encoding="utf-8")
    content = path.read_text(encoding=cfg.io.csv_encoding)
    assert content.startswith("target_chembl_id;pref_name")

    meta = path.with_name(path.name + ".meta.yaml")
    assert meta.exists()
