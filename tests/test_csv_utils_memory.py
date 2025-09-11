from __future__ import annotations

import time
from pathlib import Path
import multiprocessing as mp

import pandas as pd
import psutil

from library.csv_utils import write_csv_deterministic


def _worker(path: str, use_copy: bool, n: int) -> None:
    """Helper process writing a DataFrame to CSV."""
    df = pd.DataFrame({"a": range(n), "b": range(n)})
    if use_copy:
        df = df.copy()
    write_csv_deterministic(df, Path(path))


def _peak_memory(n: int, use_copy: bool, path: Path) -> int:
    """Return peak RSS while running ``write_csv_deterministic``."""
    proc = mp.Process(target=_worker, args=(str(path), use_copy, n))
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


def test_write_csv_deterministic_memory_usage(tmp_path: Path) -> None:
    """Writing without an extra copy uses less memory."""
    n = 1_000_000
    peak_no_copy = _peak_memory(n, False, tmp_path / "no_copy.csv")
    peak_with_copy = _peak_memory(n, True, tmp_path / "with_copy.csv")
    # ``peak_with_copy`` should not be smaller as it includes an additional
    # DataFrame allocation.
    assert peak_with_copy >= peak_no_copy
