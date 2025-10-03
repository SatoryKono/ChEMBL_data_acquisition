"""Regression tests for streamed CSV exports."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import pytest


from library.common.csv_utils import write_csv_deterministic
from library.utils.cli_tools import csv_utils_main as cli


def test_write_csv_deterministic_respects_chunksize(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Chunked exports forward the configured size to :meth:`DataFrame.to_csv`."""

    df = pd.DataFrame({"id": [2, 1], "value": ["b", "a"]})
    recorded: dict[str, object] = {}
    original_to_csv = pd.DataFrame.to_csv

    def capture_chunksize(self, *args, **kwargs):
        recorded["chunksize"] = kwargs.get("chunksize")
        return original_to_csv(self, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_csv", capture_chunksize)

    path = tmp_path / "chunked.csv"
    write_csv_deterministic(df, path, key_cols=["id"], chunksize=1)

    assert path.exists()
    assert recorded["chunksize"] == 1


def test_write_csv_deterministic_uses_chunksize(tmp_path: Path, monkeypatch) -> None:
    df = pd.DataFrame({"id": [3, 1, 2], "flag": [True, False, True]})

    recorded: dict[str, object] = {}
    original = pd.DataFrame.to_csv

    def spy(self, *args, **kwargs):  # type: ignore[override]
        recorded["chunksize"] = kwargs.get("chunksize")
        return original(self, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_csv", spy)

    out_path = tmp_path / "stream.csv"
    write_csv_deterministic(
        df,
        out_path,
        col_order=["id", "flag"],
        key_cols=["id"],
        chunksize=64,
    )

    assert recorded["chunksize"] == 64
    assert out_path.exists()


def test_cli_preserves_leading_zeroes(tmp_path: Path) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id,value\n001,a\n010,b\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"

    rc = cli.main([
        "--input",
        str(input_csv),
        "--final-out",
        str(output_csv),
        "--key-cols",
        "id",
    ])

    assert rc == 0
    contents = output_csv.read_text(encoding="utf8").splitlines()
    assert contents[1].startswith("001,")
    assert contents[2].startswith("010,")
