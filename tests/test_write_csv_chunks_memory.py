from __future__ import annotations

import multiprocessing as mp
import time
from collections.abc import Iterable
from pathlib import Path

import pandas as pd
import pytest

from library.csv_utils import write_csv_chunks_deterministic

psutil = pytest.importorskip("psutil")


def _worker(path: str, rows: int, chunks: int) -> None:
    df = pd.DataFrame({"a": range(rows), "b": range(rows)})

    def gen() -> Iterable[pd.DataFrame]:
        for _ in range(chunks):
            yield df.copy()

    write_csv_chunks_deterministic(
        gen(),
        Path(path),
        key_cols=["a"],
        merge_chunksize=50_000,
    )


def _peak_memory(path: Path, rows: int, chunks: int) -> int:
    proc = mp.Process(target=_worker, args=(str(path), rows, chunks))
    proc.start()
    ps_proc = psutil.Process(proc.pid)
    peak = 0
    while proc.is_alive():
        try:
            mem = ps_proc.memory_info().rss
            if mem > peak:
                peak = mem
        except psutil.NoSuchProcess:
            break
        time.sleep(0.05)
    proc.join()
    return peak


def test_write_csv_chunks_memory_usage(tmp_path: Path) -> None:
    rows = 200_000
    peak_single = _peak_memory(tmp_path / "one.csv", rows, 1)
    peak_many = _peak_memory(tmp_path / "many.csv", rows, 5)
    assert peak_many <= peak_single * 1.5
