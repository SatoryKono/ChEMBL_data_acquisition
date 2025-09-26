"""Regression tests for streamed CSV exports."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.csv_utils import write_csv_deterministic


def test_write_csv_deterministic_uses_chunksize(
    tmp_path: Path, monkeypatch
) -> None:
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
