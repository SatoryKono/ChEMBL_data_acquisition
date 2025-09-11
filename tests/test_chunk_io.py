"""Tests for chunked CSV processing with checkpoints."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.chunk_io import process_csv_chunks
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
