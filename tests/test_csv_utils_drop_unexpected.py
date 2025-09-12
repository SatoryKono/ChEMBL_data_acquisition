from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library import csv_utils


def test_drop_unexpected_columns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    df = pd.DataFrame({"a": [1], "b": [2], "c": [3]})
    out = tmp_path / "out.csv"
    messages: list[str] = []
    monkeypatch.setattr(
        csv_utils.logger,
        "warning",
        lambda msg, *args, **kwargs: messages.append(msg),
    )
    csv_utils.write_csv_deterministic(
        df,
        out,
        col_order=["a", "b"],
        key_cols=["a"],
        drop_unexpected_cols=True,
    )
    result = pd.read_csv(out)
    assert list(result.columns) == ["a", "b"]
    assert "unexpected_columns_dropped" in messages
