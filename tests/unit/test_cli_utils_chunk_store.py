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

    paths = list(store._paths)  # type: ignore[attr-defined]
    assert len(paths) == 1
    assert paths[0].suffix == ".pkl"

    store.cleanup()


@pytest.mark.unit
def test_parquet_chunk_store__append_tracks_single_frame_row_count():
    """Appending a single chunk persists it once and updates counters."""

    store = cli_utils._ParquetChunkStore()

    frame = pd.DataFrame({"col": [1, 2, 3]})
    store.append(frame)

    frames = list(store.iter_frames())
    assert len(frames) == 1
    pd.testing.assert_frame_equal(frames[0], frame)

    assert store.row_count == len(frame)

    store.cleanup()

