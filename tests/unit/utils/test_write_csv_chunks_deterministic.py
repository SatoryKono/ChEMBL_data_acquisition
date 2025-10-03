from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.common.csv_utils import write_csv_chunks_deterministic, write_csv_deterministic


def test_write_csv_chunks_deterministic(tmp_path: Path) -> None:
    """Chunk iterator writes must match full DataFrame output."""
    rows = 5000
    df = pd.DataFrame(
        {
            "a": range(rows),
            "b": [i % 2 == 0 for i in range(rows)],
            "d": pd.date_range("2020-01-01", periods=rows, freq="D"),
            "f": [i / 3 for i in range(rows)],
            "chembl_id": pd.Series((f"{i:05d}" for i in range(rows)), dtype="string"),
        }
    )
    shuffled = df.sample(frac=1, random_state=1).reset_index(drop=True)
    path_full = tmp_path / "full.csv"
    path_chunks = tmp_path / "chunks.csv"
    write_csv_deterministic(
        shuffled,
        path_full,
        col_order=["chembl_id", "a", "b", "d", "f"],
        key_cols=["a"],
    )
    chunk_iter = (shuffled.iloc[i : i + 1000] for i in range(0, len(shuffled), 1000))
    write_csv_chunks_deterministic(
        chunk_iter,
        path_chunks,
        col_order=["chembl_id", "a", "b", "d", "f"],
        key_cols=["a"],
        chunksize=1000,
    )
    assert path_full.read_bytes() == path_chunks.read_bytes()


def test_write_csv_chunks_deterministic_missing_keys(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "name": ["beta", "alpha", "gamma", "delta"],
            "value": [2, 1, 3, 1],
        }
    )
    shuffled = df.sample(frac=1, random_state=5).reset_index(drop=True)
    path_full = tmp_path / "full_fallback.csv"
    path_chunks = tmp_path / "chunks_fallback.csv"

    write_csv_deterministic(
        shuffled,
        path_full,
        col_order=["value", "name"],
        key_cols=["chembl_id"],
    )

    chunk_iter = (shuffled.iloc[i : i + 2] for i in range(0, len(shuffled), 2))
    write_csv_chunks_deterministic(
        chunk_iter,
        path_chunks,
        col_order=["value", "name"],
        key_cols=["chembl_id"],
        chunksize=2,
        merge_chunksize=1,
    )

    assert path_full.read_bytes() == path_chunks.read_bytes()


def test_write_csv_chunks_preserves_leading_zeros(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "chembl_id": pd.Series(["0001", "0100", "1000"], dtype="string"),
            "value": [3, 2, 1],
        }
    )

    shuffled = df.sample(frac=1, random_state=7).reset_index(drop=True)
    path_full = tmp_path / "full_leading.csv"
    path_chunks = tmp_path / "chunks_leading.csv"

    write_csv_deterministic(
        shuffled,
        path_full,
        col_order=["chembl_id", "value"],
        key_cols=["chembl_id"],
    )

    chunk_iter = (shuffled.iloc[i : i + 1] for i in range(0, len(shuffled)))
    write_csv_chunks_deterministic(
        chunk_iter,
        path_chunks,
        col_order=["chembl_id", "value"],
        key_cols=["chembl_id"],
        chunksize=1,
        merge_chunksize=1,
    )

    assert path_full.read_bytes() == path_chunks.read_bytes()


def test_write_csv_chunks_respects_file_window(tmp_path: Path) -> None:
    rows = 60
    df = pd.DataFrame(
        {
            "chembl_id": pd.Series(
                (f"CHEMBL{i:05d}" for i in range(rows)), dtype="string"
            ),
            "value": list(range(rows)),
        }
    )
    shuffled = df.sample(frac=1, random_state=11).reset_index(drop=True)
    path_full = tmp_path / "full_window.csv"
    path_chunks = tmp_path / "chunks_window.csv"

    write_csv_deterministic(
        shuffled,
        path_full,
        col_order=["chembl_id", "value"],
        key_cols=["chembl_id"],
    )

    chunk_iter = (shuffled.iloc[i : i + 10] for i in range(0, len(shuffled), 10))
    write_csv_chunks_deterministic(
        chunk_iter,
        path_chunks,
        col_order=["chembl_id", "value"],
        key_cols=["chembl_id"],
        chunksize=10,
        merge_chunksize=5,
        merge_file_window=2,
    )

    assert path_full.read_bytes() == path_chunks.read_bytes()
