from __future__ import annotations

import hashlib
import logging
import subprocess
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from hypothesis.extra.pandas import column, data_frames, range_indexes

from library.config import Config
from library.csv_utils import _git_sha, sha256_file, write_csv_deterministic


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
        df, path, col_order=["a", "b", "d", "f"], key_cols=["a"], cfg=Config()
    )
    assert result == path
    text = path.read_text(encoding="utf-8-sig")
    assert text == ("a,b,d,f\n1,false,2020-01-01,2.34568\n2,true,2020-01-02,1.23457\n")
    assert Path(str(path) + ".meta.yaml").exists()


def test_write_csv_deterministic_hash(tmp_path: Path) -> None:
    """Generated CSV has stable SHA-256 hash."""

    df = pd.DataFrame(
        {
            "b": [True, False],
            "a": [2, 1],
            "d": [pd.Timestamp("2020-01-02"), pd.Timestamp("2020-01-01")],
            "f": [1.23456789, 2.3456789],
        }
    )
    path = tmp_path / "out.csv"
    write_csv_deterministic(df, path, col_order=["a", "b", "d", "f"], key_cols=["a"])

    expected_bytes = (
        "a,b,d,f\n1,false,2020-01-01,2.34568\n2,true,2020-01-02,1.23457\n"
    ).encode("utf-8-sig")
    expected_hash = hashlib.sha256(expected_bytes).hexdigest()

    assert sha256_file(path) == expected_hash


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


@settings(deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(
    df=data_frames(
        columns=[
            column("a", elements=st.integers(min_value=0, max_value=10), unique=True),
            column("b", elements=st.booleans()),
            column(
                "d",
                elements=st.datetimes(
                    min_value=pd.Timestamp("1970-01-01").to_pydatetime(),
                    max_value=pd.Timestamp("2100-12-31").to_pydatetime(),
                ),
            ),
            column(
                "f",
                elements=st.floats(allow_nan=False, allow_infinity=False, width=32),
            ),
        ],
        index=range_indexes(min_size=1, max_size=5),
    )
)
def test_write_csv_deterministic_hypothesis(tmp_path: Path, df: pd.DataFrame) -> None:
    """Repeated writes yield identical outputs for random data."""

    path1 = tmp_path / "first.csv"
    path2 = tmp_path / "second.csv"
    df1 = df.sample(frac=1).reset_index(drop=True)
    df2 = df.sample(frac=1).reset_index(drop=True)
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


def test_git_sha_timeout_returns_unknown_and_logs_warning(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """_git_sha returns 'unknown' and logs a warning on timeout."""

    with patch(
        "library.csv_utils.subprocess.run",
        side_effect=subprocess.TimeoutExpired(
            cmd=["git", "rev-parse", "HEAD"], timeout=5
        ),
    ):
        with caplog.at_level(logging.WARNING):
            result = _git_sha()
    assert result == "unknown"
    assert "timed out" in caplog.text
