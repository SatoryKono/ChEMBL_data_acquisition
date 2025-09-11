from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.csv_utils import write_csv_deterministic


def test_drop_unexpected_columns(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    df = pd.DataFrame({"a": [1], "b": [2], "c": [3]})
    out = tmp_path / "out.csv"
    with caplog.at_level("WARNING"):
        write_csv_deterministic(
            df,
            out,
            col_order=["a", "b"],
            drop_unexpected_cols=True,
        )
    result = pd.read_csv(out)
    assert list(result.columns) == ["a", "b"]
    assert "Dropping unexpected columns" in caplog.text
