"""Tests for chunked CSV processing with checkpoints."""

from __future__ import annotations

import gc
import warnings
from pathlib import Path

import pandas as pd

from library.io.chunk_io import process_csv_chunks, read_csv_chunks
from library.config import IoCfg


def test_process_csv_chunks_resume(tmp_path: Path) -> None:
    """Processing with a checkpoint should resume from last offset."""
    cfg = IoCfg()
    input_path = tmp_path / "input.csv"
    df = pd.DataFrame({"a": range(200)})
    df.to_csv(input_path, index=False)

    output_path = tmp_path / "output.csv"
    checkpoint = tmp_path / "checkpoint.json"

    # First run processes only 100 rows
    rows_written = process_csv_chunks(
        input_path,
        output_path,
        cfg=cfg,
        chunk_size=50,
        limit=100,
        checkpoint_path=checkpoint,
    )
    assert rows_written == 100
    interim = pd.read_csv(output_path)
    assert len(interim) == 100

    # Second run resumes and completes the remaining rows
    rows_written = process_csv_chunks(
        input_path,
        output_path,
        cfg=cfg,
        chunk_size=50,
        checkpoint_path=checkpoint,
    )
    assert rows_written == 200
    result = pd.read_csv(output_path)
    assert result.equals(df)


def test_read_csv_chunks_no_resource_warning(tmp_path: Path) -> None:
    """``read_csv_chunks`` should not emit ``ResourceWarning``."""

    cfg = IoCfg()
    path = tmp_path / "input.csv"
    pd.DataFrame({"a": range(10)}).to_csv(path, index=False)

    with warnings.catch_warnings(record=True) as warns:
        warnings.simplefilter("always", ResourceWarning)
        gen = read_csv_chunks(path, cfg=cfg, chunk_size=5)
        next(gen)
        del gen
        gc.collect()

    assert not any(w.category is ResourceWarning for w in warns)
