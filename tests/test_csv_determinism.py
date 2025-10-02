from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.common.csv_utils import sha256_file, write_csv_deterministic


def test_write_csv_deterministic_hash(tmp_path: Path) -> None:
    """Writing the same DataFrame twice yields identical SHA-256 hashes."""
    rows = 100
    df = pd.DataFrame(
        {
            "a": range(rows),
            "b": [i % 2 == 0 for i in range(rows)],
            "d": pd.date_range("2020-01-01", periods=rows, freq="D"),
            "f": [i / 3 for i in range(rows)],
        }
    )
    shuffled = df.sample(frac=1, random_state=42).reset_index(drop=True)
    first_path = tmp_path / "first.csv"
    second_path = tmp_path / "second.csv"
    write_csv_deterministic(
        shuffled.copy(),
        first_path,
        col_order=["a", "b", "d", "f"],
        key_cols=["a"],
    )
    write_csv_deterministic(
        shuffled.copy(),
        second_path,
        col_order=["a", "b", "d", "f"],
        key_cols=["a"],
    )
    first_hash = sha256_file(first_path)
    second_hash = sha256_file(second_path)
    assert first_hash == second_hash
