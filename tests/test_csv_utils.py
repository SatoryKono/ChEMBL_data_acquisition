from __future__ import annotations

import hashlib
import string
import subprocess
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
import pytest

pytest.importorskip("hypothesis")
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from hypothesis.extra.pandas import column, data_frames, range_indexes

import library.git_utils as git_utils
from library.config import Config
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
    cfg = Config(api={"user_agent": "test@example.org"})
    result = write_csv_deterministic(
        df, path, col_order=["a", "b", "d", "f"], key_cols=["a"], cfg=cfg
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


def test_write_csv_deterministic_duplicate_columns(tmp_path: Path) -> None:
    """Duplicate column labels raise a ValueError."""

    df = pd.DataFrame([[1, 2]], columns=["a", "a"])
    path = tmp_path / "dup.csv"
    msg = r"Duplicate column names found: \['a'\]"
    with pytest.raises(ValueError, match=msg):
        write_csv_deterministic(df, path, key_cols=["a"])


def test_default_sorting_and_order(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "b": [False, True],
            "a": ["y", "x"],
        }
    )
    path = tmp_path / "out.csv"
    write_csv_deterministic(df, path, key_cols=["a"])
    assert path.read_text(encoding="utf-8-sig") == "a,b\nx,true\ny,false\n"


def test_write_csv_deterministic_empty_dataframe_golden(tmp_path: Path) -> None:
    """Empty frames retain deterministic header order."""

    df = pd.DataFrame(columns=["b", "a"])
    path = tmp_path / "empty.csv"
    write_csv_deterministic(df, path, key_cols=sorted(df.columns))

    expected = (
        Path(__file__).parent
        / "data"
        / "golden"
        / "empty_with_header.csv"
    ).read_text(encoding="utf-8-sig")
    assert path.read_text(encoding="utf-8-sig") == expected


def test_write_csv_deterministic_empty_no_columns(tmp_path: Path) -> None:
    """DataFrames without columns can still be written deterministically."""

    df = pd.DataFrame()
    path = tmp_path / "empty.csv"
    write_csv_deterministic(df, path, key_cols=[])

    expected = (
        Path(__file__).parent
        / "data"
        / "golden"
        / "empty_no_columns.csv"
    ).read_text(encoding="utf-8-sig")
    assert path.read_text(encoding="utf-8-sig") == expected


def test_write_csv_deterministic_requires_key_columns(tmp_path: Path) -> None:
    df = pd.DataFrame({"a": [1], "b": [2]})
    path = tmp_path / "out.csv"
    msg = "key_cols must contain at least one column"
    with pytest.raises(ValueError, match=msg):
        write_csv_deterministic(df, path, key_cols=[])


def test_write_csv_deterministic_none_key_cols(tmp_path: Path) -> None:
    df = pd.DataFrame({"a": [1], "b": [2]})
    path = tmp_path / "out.csv"
    with pytest.raises(ValueError, match="key_cols must be provided"):
        write_csv_deterministic(df, path, key_cols=None)


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


def test_write_csv_deterministic_hash_stable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Repeated calls maintain identical data and metadata hashes."""

    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    path = tmp_path / "out.csv"

    # Freeze timestamp to keep metadata deterministic across runs
    fixed_now = datetime(2024, 1, 1)
    monkeypatch.setattr("library.io.datetime", SimpleNamespace(now=lambda: fixed_now))

    write_csv_deterministic(df.copy(), path, key_cols=sorted(df.columns))
    first_hash = sha256_file(path)
    first_meta_hash = sha256_file(Path(str(path) + ".meta.yaml"))

    write_csv_deterministic(df.copy(), path, key_cols=sorted(df.columns))
    second_hash = sha256_file(path)
    second_meta_hash = sha256_file(Path(str(path) + ".meta.yaml"))

    assert first_hash == second_hash
    assert first_meta_hash == second_meta_hash


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


@settings(deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(data=st.data())
def test_write_csv_deterministic_random_types(
    tmp_path: Path, data: st.DataObject
) -> None:
    """DataFrames with mixed dtypes yield deterministic CSV output."""

    element_strategies = [
        st.integers(min_value=-10, max_value=10),
        st.floats(allow_nan=False, allow_infinity=False, width=32),
        st.booleans(),
        st.text(alphabet=string.ascii_letters, min_size=1, max_size=5),
        st.datetimes(
            min_value=pd.Timestamp("1970-01-01").to_pydatetime(),
            max_value=pd.Timestamp("2100-12-31").to_pydatetime(),
        ),
    ]

    columns = data.draw(
        st.lists(
            st.builds(
                column,
                name=st.text(
                    alphabet=string.ascii_letters,
                    min_size=1,
                    max_size=5,
                ).filter(lambda n: n != "a"),
                elements=st.one_of(*[st.just(es) for es in element_strategies]),
            ),
            min_size=1,
            max_size=5,
            unique_by=lambda c: c.name,
        )
    )
    df = data.draw(
        data_frames(columns=columns, index=range_indexes(min_size=1, max_size=5))
    )
    df.insert(0, "a", range(len(df)))

    path1 = tmp_path / "first.csv"
    path2 = tmp_path / "second.csv"
    df1 = df.sample(frac=1).reset_index(drop=True)
    df2 = df.sample(frac=1).reset_index(drop=True)
    write_csv_deterministic(df1.copy(), path1, key_cols=["a"])
    write_csv_deterministic(df2.copy(), path2, key_cols=["a"])
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
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """_git_sha returns 'unknown' and logs a warning on timeout."""

    git_utils._git_sha.cache_clear()
    monkeypatch.setattr(git_utils, "_read_head_sha", lambda *_: None)
    with patch(
        "library.git_utils.subprocess.run",
        side_effect=subprocess.CalledProcessError(
            returncode=1,
            cmd=["git", "rev-parse", "HEAD"],
            stderr="fatal: simulated error",
        ),
    ):
        records: list[tuple[str, dict[str, str] | None]] = []

        def fake_warning(
            event: str,
            *args: object,
            extra: dict[str, str] | None = None,
            **kwargs: object,
        ) -> None:
            records.append((event, extra))

        monkeypatch.setattr(git_utils.logger, "warning", fake_warning)
        result = git_utils._git_sha()

    assert result == "UNKNOWN"
    assert any(
        event == "git_sha_unavailable" and extra is not None and "error" in extra
        for event, extra in records
    )
