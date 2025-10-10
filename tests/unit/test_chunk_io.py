"""Unit tests for :mod:`library.common.chunk_io`."""

from __future__ import annotations

from library.common.chunk_io import process_csv_chunks
from library.config import IoCfg


def test_process_csv_chunks__writes_unix_newlines(tmp_path):
    """Ensure chunked CSV copy always emits Unix line endings."""
    input_path = tmp_path / "input.csv"
    with input_path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("col1,col2\r\n1,2\r\n3,4\r\n")

    output_path = tmp_path / "output.csv"
    cfg = IoCfg()

    rows = process_csv_chunks(
        input_path,
        output_path,
        cfg=cfg,
        chunk_size=1,
        ensure_directory=False,
    )

    assert rows == 2
    content = output_path.read_bytes()
    assert content == b"col1,col2\n1,2\n3,4\n"
