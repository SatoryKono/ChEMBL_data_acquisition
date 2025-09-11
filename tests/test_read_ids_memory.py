from __future__ import annotations

from collections import deque
from pathlib import Path
import tracemalloc

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
    deque(ids, maxlen=0)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak < 10 * 1024 * 1024
