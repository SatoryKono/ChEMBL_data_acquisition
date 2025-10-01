from __future__ import annotations

import multiprocessing as mp
import time
from pathlib import Path

import pandas as pd
import pytest

from library.config import Config
from library.csv_utils import write_csv_deterministic

psutil = pytest.importorskip("psutil")


def _worker(path: str, use_copy: bool, n: int) -> None:
    """Helper process writing a DataFrame to CSV."""
    df = pd.DataFrame({"a": range(n), "b": range(n)})
    if use_copy:
        df = df.copy()
    cfg = Config()
    cfg.io.exist_ok = False
    write_csv_deterministic(df, Path(path), key_cols=["a", "b"], cfg=cfg)


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
    # ``peak_with_copy`` should not be significantly smaller as it includes an
    # additional DataFrame allocation. Allow a small tolerance for measurement
    # noise introduced by RSS sampling.
    assert peak_with_copy >= 0.99 * peak_no_copy


def test_write_csv_deterministic_memory_usage_large(tmp_path: Path) -> None:
    """Handles DataFrames larger than one million rows."""
    n = 1_200_000
    peak_no_copy = _peak_memory(n, False, tmp_path / "large_no_copy.csv")
    peak_with_copy = _peak_memory(n, True, tmp_path / "large_with_copy.csv")
    # Allow minor fluctuations in memory measurements on different platforms
    assert peak_with_copy >= 0.9 * peak_no_copy
