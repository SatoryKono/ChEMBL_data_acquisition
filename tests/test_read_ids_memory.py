from __future__ import annotations

import tracemalloc
from collections import deque
from collections.abc import Iterator
from pathlib import Path

import pytest

from library import io
from library.config import IoCfg


def test_read_ids_streams_without_memory_growth(tmp_path: Path) -> None:
    path = tmp_path / "ids.csv"
    with path.open("w", newline="") as fh:
        fh.write("id\n")
        for i in range(100000):
            fh.write(f"{i}\n")

    tracemalloc.start()
    ids = io.read_ids(path, column="id", cfg=IoCfg())
    assert isinstance(ids, Iterator)
    assert iter(ids) is ids
    deque(ids, maxlen=0)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak < 10 * 1024 * 1024


def test_read_ids_logs_dropped_markers_on_close(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "ids.csv"
    with path.open("w", newline="") as fh:
        fh.write("id\n")
        fh.write("CHEMBL1\n")
        fh.write("NA\n")

    calls: list[tuple[str, tuple[object, ...], dict[str, object]]] = []

    from library.io import readers as io_readers

    def _record(event: str, *args: object, **data: object) -> None:
        calls.append((event, args, data))

    monkeypatch.setattr(io_readers.logger, "warning", _record)
    ids = io.read_ids(path, column="id", cfg=IoCfg(), na_markers=["NA"], keep_na_markers=False)
    assert next(ids) == "CHEMBL1"
    with pytest.raises(StopIteration):
        next(ids)

    assert (
        "read_ids_dropped_na_markers",
        (),
        {
            "path": str(path),
            "column": "id",
            "dropped_total": 1,
            "dropped_ids": ["NA"],
        },
    ) in calls
