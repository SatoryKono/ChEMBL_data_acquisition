from __future__ import annotations

from pathlib import Path
import hashlib

import pandas as pd
import pytest

from library.csv_utils import sha256_file, write_csv_deterministic


def test_write_csv_deterministic(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "b": [True, False],
            "a": [2, 1],
            "d": [pd.Timestamp("2020-01-02"), pd.Timestamp("2020-01-01")],
            "f": [1.23456789, 2.3456789],
        }
    )
    path = tmp_path / "out.csv"
    result = write_csv_deterministic(
        df, path, col_order=["a", "b", "d", "f"], key_cols=["a"]
    )
    assert result == path
    text = path.read_text(encoding="utf-8-sig")
    assert text == (
        "a,b,d,f\n" "1,false,2020-01-01,2.34568\n" "2,true,2020-01-02,1.23457\n"
    )


def test_default_sorting_and_order(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "b": [False, True],
            "a": ["y", "x"],
        }
    )
    path = tmp_path / "out.csv"
    write_csv_deterministic(df, path)
    assert path.read_text(encoding="utf-8-sig") == "a,b\nx,true\ny,false\n"


def test_deterministic_writes_identical_bytes(tmp_path: Path) -> None:
    """Ensure deterministic writes produce identical files."""

    data_path = Path(__file__).parent / "data" / "csv_utils_input.csv"
    df = pd.read_csv(data_path, parse_dates=["d"])
    df1 = df.sample(frac=1, random_state=1).reset_index(drop=True)
    df2 = df.sample(frac=1, random_state=2).reset_index(drop=True)

    path1 = tmp_path / "first.csv"
    path2 = tmp_path / "second.csv"

    write_csv_deterministic(df1, path1, col_order=["a", "b", "d", "f"], key_cols=["a"])
    write_csv_deterministic(df2, path2, col_order=["a", "b", "d", "f"], key_cols=["a"])

    assert path1.read_bytes() == path2.read_bytes()


def test_sha256_file(tmp_path: Path) -> None:
    data = b"hello world"
    file_path = tmp_path / "data.txt"
    file_path.write_bytes(data)
    expected = hashlib.sha256(data).hexdigest()
    assert sha256_file(file_path) == expected


def test_sha256_file_missing(tmp_path: Path) -> None:
    missing = tmp_path / "missing.csv"
    with pytest.raises(FileNotFoundError, match=str(missing)):
        sha256_file(missing)
