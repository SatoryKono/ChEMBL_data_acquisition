from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.csv_utils import write_csv_deterministic


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
