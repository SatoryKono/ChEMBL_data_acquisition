"""Chunked CSV processing utilities with checkpoint support.

This module provides helper functions to read and write CSV files in
manageable chunks. A simple checkpointing mechanism records the number of
processed rows enabling pipelines to resume from the last successful chunk.
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from pathlib import Path
from typing import Literal

import pandas as pd

from ..config import IoCfg
from .log import logger


def _read_checkpoint(path: Path) -> int:
    """Return number of rows already processed.

    Parameters
    ----------
    path:
        Location of the checkpoint file.

    Returns
    -------
    int
        Number of rows previously written. ``0`` if the checkpoint does not
        exist or is malformed.
    """
    try:
        with path.open("r", encoding="utf8") as fh:
            data = json.load(fh)
        rows = int(data.get("rows", 0))
        if rows < 0:
            raise ValueError
        return rows
    except FileNotFoundError:
        return 0
    except (ValueError, json.JSONDecodeError) as exc:
        logger.warning("invalid_checkpoint", path=str(path), error=str(exc))
        return 0


def _write_checkpoint(path: Path, rows: int) -> None:
    """Persist ``rows`` to ``path`` as JSON."""
    tmp = path.with_suffix(".tmp")
    with tmp.open("w", encoding="utf8") as fh:
        json.dump({"rows": rows}, fh)
    tmp.replace(path)


def read_csv_chunks(
    path: str | Path,
    *,
    cfg: IoCfg,
    chunk_size: int,
    limit: int | None = None,
    skiprows: int = 0,
    sep: str | None = None,
    encoding: str | None = None,
) -> Iterator[pd.DataFrame]:
    """Yield chunks from ``path``.

    Parameters
    ----------
    path:
        CSV input file.
    cfg:
        I/O defaults for CSV handling.
    chunk_size:
        Maximum number of rows per yielded chunk. Must be positive.
    limit:
        Optional global row limit across all chunks.
    skiprows:
        Number of data rows to skip from the start (for resume).
    sep:
        Field delimiter overriding ``cfg.csv_sep``.
    encoding:
        Character encoding overriding ``cfg.csv_encoding``.

    Yields
    ------
    pandas.DataFrame
        Consecutive chunks of the input table.
    """
    if chunk_size <= 0:
        msg = f"chunk_size must be positive, got {chunk_size}"
        raise ValueError(msg)

    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding

    # ``skiprows`` expects an iterable of row indices starting at 0 including
    # the header. We skip the header plus ``skiprows`` data rows.
    to_skip = range(1, skiprows + 1) if skiprows else None
    reader = pd.read_csv(
        path,
        sep=sep,
        encoding=encoding,
        chunksize=chunk_size,
        skiprows=to_skip,
    )
    emitted = 0
    try:
        for chunk in reader:
            if limit is not None and emitted >= limit:
                break
            if limit is not None:
                remaining = limit - emitted
                if len(chunk) > remaining:
                    yield chunk.iloc[:remaining]
                    break
            emitted += len(chunk)
            yield chunk
    finally:
        reader.close()


def process_csv_chunks(
    input_path: str | Path,
    output_path: str | Path,
    *,
    cfg: IoCfg,
    chunk_size: int,
    limit: int | None = None,
    checkpoint_path: Path | None = None,
    sep: str | None = None,
    encoding: str | None = None,
    ensure_directory: bool = False,
    line_terminator: str = "\n",
) -> int:
    """Copy ``input_path`` to ``output_path`` using chunked I/O.

    Existing data and checkpoints allow resuming interrupted runs. The number of
    rows written is returned.

    Parameters
    ----------
    input_path:
        Source CSV file.
    output_path:
        Destination CSV file.
    cfg:
        I/O configuration providing defaults.
    chunk_size:
        Number of rows processed per chunk.
    limit:
        Optional maximum number of rows to process in total. ``None`` processes
        the entire file.
    checkpoint_path:
        Location of a checkpoint file storing processed row counts.
    sep, encoding:
        CSV formatting options overriding ``cfg`` defaults.
    ensure_directory:
        When ``True``, create the parent directory of ``output_path`` if it does
        not exist and :attr:`IoCfg.exist_ok` permits directory creation. When
        ``False``, callers must ensure the directory is available.
    line_terminator:
        Line ending to use when serialising CSV data. Defaults to a Unix newline
        to guarantee consistent output across platforms.

    Returns
    -------
    int
        Total number of rows written after completion.
    """
    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding

    output_path = Path(output_path)
    if ensure_directory:
        parent = output_path.parent
        if parent.exists():
            if not parent.is_dir():
                msg = f"{parent} is not a directory"
                raise NotADirectoryError(msg)
        else:
            if cfg.exist_ok:
                parent.mkdir(parents=True, exist_ok=True)
            else:
                msg = f"{parent} does not exist"
                raise FileNotFoundError(msg)

    processed = 0
    header_written = False
    if checkpoint_path is not None:
        processed = _read_checkpoint(checkpoint_path)
        header_written = output_path.exists() and processed > 0
        if processed:
            logger.info("resume_from_row", row=processed)

    rows_written = processed
    for chunk in read_csv_chunks(
        input_path,
        cfg=cfg,
        chunk_size=chunk_size,
        limit=None if limit is None else max(limit - rows_written, 0),
        skiprows=processed,
        sep=sep,
        encoding=encoding,
    ):
        mode: Literal["a", "w"] = "a" if header_written else "w"
        chunk.to_csv(
            path_or_buf=output_path,
            mode=mode,
            index=False,
            header=not header_written,
            sep=sep,
            encoding=encoding,
            lineterminator=line_terminator,
        )
        header_written = True
        rows_written += len(chunk)
        if checkpoint_path is not None:
            _write_checkpoint(checkpoint_path, rows_written)
        if limit is not None and rows_written >= limit:
            break

    return rows_written
