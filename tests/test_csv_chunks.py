"""Tests for chunked deterministic CSV writing."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.csv_utils import write_csv_deterministic


def test_write_csv_deterministic_chunks(tmp_path: Path) -> None:
    """Chunked and unchunked writes must yield identical files."""

    rows = 5000
    df = pd.DataFrame(
        {
            "a": range(rows),
            "b": [i % 2 == 0 for i in range(rows)],
            "d": pd.date_range("2020-01-01", periods=rows, freq="D"),
            "f": [i / 3 for i in range(rows)],
        }
    )
    shuffled = df.sample(frac=1, random_state=1).reset_index(drop=True)
    path1 = tmp_path / "unchunked.csv"
    path2 = tmp_path / "chunked.csv"
    write_csv_deterministic(
        shuffled,
        path1,
        col_order=["a", "b", "d", "f"],
        key_cols=["a"],
    )
    write_csv_deterministic(
        shuffled,
        path2,
        col_order=["a", "b", "d", "f"],
        key_cols=["a"],
        chunksize=1000,
    )
    assert path1.read_bytes() == path2.read_bytes()
