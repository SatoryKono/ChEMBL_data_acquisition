from __future__ import annotations

import importlib

import pandas as pd
import pytest

from library import cli_utils


@pytest.mark.unit
def test_parquet_chunk_store__falls_back_to_pickle_when_engine_missing(monkeypatch):
    """Ensure parquet chunk store degrades gracefully without optional deps."""

    real_import_module = importlib.import_module

    def fake_import_module(name: str, package: str | None = None):
        if name in {"pyarrow", "fastparquet"}:
            raise ImportError(name)
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", fake_import_module)

    real_to_parquet = pd.DataFrame.to_parquet
    real_to_pickle = pd.DataFrame.to_pickle
    pickle_calls: list[str] = []

    def fail_to_parquet(self, *args, **kwargs):  # type: ignore[override]
        raise AssertionError("to_parquet should not be invoked when unavailable")

    def record_to_pickle(self, path, *args, **kwargs):  # type: ignore[override]
        pickle_calls.append(str(path))
        return real_to_pickle(self, path, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_parquet", fail_to_parquet)
    monkeypatch.setattr(pd.DataFrame, "to_pickle", record_to_pickle)

    store = cli_utils._ParquetChunkStore()

    frame = pd.DataFrame({"col": [1, 2, 3]})
    store.append(frame)

    assert len(pickle_calls) == 1

    loaded_frames = list(store.iter_frames())
    assert len(loaded_frames) == 1
    pd.testing.assert_frame_equal(loaded_frames[0], frame)

    chunks = list(store._chunks)  # type: ignore[attr-defined]
    assert len(chunks) == 1
    assert chunks[0].path.suffix == ".pkl"
    assert chunks[0].kind == "pickle"

    store.cleanup()


@pytest.mark.unit
def test_parquet_chunk_store__reads_mixed_backends(monkeypatch):
    """First chunk stored as parquet, second falls back to pickle."""

    real_to_pickle = pd.DataFrame.to_pickle
    real_read_pickle = pd.read_pickle

    parquet_writes: list[str] = []
    parquet_reads: list[str] = []
    pickle_reads: list[str] = []
    call_count = {"to_parquet": 0}

    def fake_to_parquet(self, path, *args, **kwargs):  # type: ignore[override]
        if call_count["to_parquet"] == 0:
            call_count["to_parquet"] += 1
            parquet_writes.append(str(path))
            return real_to_pickle(self, path)
        raise ImportError("parquet engine unavailable")

    def fake_read_parquet(path, *args, **kwargs):  # type: ignore[override]
        parquet_reads.append(str(path))
        return real_read_pickle(path, *args, **kwargs)

    def record_read_pickle(path, *args, **kwargs):  # type: ignore[override]
        pickle_reads.append(str(path))
        return real_read_pickle(path, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_parquet", fake_to_parquet)
    monkeypatch.setattr(pd, "read_parquet", fake_read_parquet)
    monkeypatch.setattr(pd, "read_pickle", record_read_pickle)

    store = cli_utils._ParquetChunkStore()

    frame1 = pd.DataFrame({"col": [1, 2, 3]})
    frame2 = pd.DataFrame({"col": [4, 5, 6]})

    store.append(frame1)
    store.append(frame2)

    chunks = list(store._chunks)  # type: ignore[attr-defined]
    assert len(chunks) == 2
    assert chunks[0].kind == "parquet"
    assert chunks[1].kind == "pickle"

    loaded_frames = list(store.iter_frames())
    assert parquet_writes == [str(chunks[0].path)]
    assert parquet_reads == [str(chunks[0].path)]
    assert pickle_reads == [str(chunks[1].path)]
    assert len(loaded_frames) == 2
    pd.testing.assert_frame_equal(loaded_frames[0], frame1)
    pd.testing.assert_frame_equal(loaded_frames[1], frame2)

    store.cleanup()

