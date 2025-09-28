"""CSV writing helpers with deterministic output and checkpointing support."""

from __future__ import annotations

import csv
import hashlib
import heapq
import json
import os
import sys
import tempfile
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from contextlib import ExitStack
from datetime import date, datetime
from pathlib import Path
from typing import Any, Literal, TextIO, cast

import numpy as np
import pandas as pd
import yaml
from pandas.api import types as ptypes

from ..config import Config, IoCfg, _serialize_paths
from ..utils.git import git_sha
from ..utils.logging import logger
from .readers import read_csv_chunks


def _normalise_bool(series: pd.Series) -> pd.Series:
  """Return ``series`` with booleans converted to ``"true"``/``"false"``."""

  return series.map({True: "true", False: "false"}).astype("string")


def _normalise_dates(series: pd.Series) -> pd.Series:
  """Return ``series`` with dates formatted as ``YYYY-MM-DD``."""

  if ptypes.is_datetime64_any_dtype(series):
    return series.dt.strftime("%Y-%m-%d")

  if (
    ptypes.is_object_dtype(series)
    and series.dropna().map(lambda value: isinstance(value, date | datetime)).all()
  ):
    result: pd.Series = pd.to_datetime(series).dt.strftime("%Y-%m-%d")
    return result

  return series


def _resolve_sort_columns(
  df: pd.DataFrame,
  requested: Sequence[str],
  *,
  emit_warning: bool = True,
) -> list[str]:
  """Return stable row sort columns honouring ``requested`` when possible."""

  sort_cols = list(requested)
  if not sort_cols:
    if df.empty:
      return []
    raise ValueError("key_cols must contain at least one column")

  missing = [col for col in sort_cols if col not in df.columns]
  if not missing:
    return sort_cols

  available = [col for col in sort_cols if col in df.columns]
  fallback = available + [col for col in df.columns if col not in available]
  if emit_warning:
    logger.warning(
      "missing_key_columns",
      requested=sort_cols,
      missing=missing,
      fallback=fallback,
    )
  return fallback


def write_csv_deterministic(
  df: pd.DataFrame,
  path: str | Path,
  *,
  col_order: Sequence[str] | None = None,
  key_cols: Sequence[str] | None = None,
  chunksize: int | None = None,
  sort_chunksize: int | None = None,
  sep: str = ",",
  encoding: str = "utf-8-sig",
  cfg: Config | None = None,
  drop_unexpected_cols: bool = False,
) -> Path:
  """Serialise ``df`` to ``path`` as a deterministic CSV file."""

  out_path = Path(path)
  out_path.parent.mkdir(parents=True, exist_ok=True)

  if chunksize is not None and chunksize <= 0:
    raise ValueError(f"chunksize must be a positive integer, got {chunksize}")
  if sort_chunksize is not None and sort_chunksize <= 0:
    raise ValueError(f"sort_chunksize must be a positive integer, got {sort_chunksize}")

  duplicated = df.columns[df.columns.duplicated()]
  if not duplicated.empty:
    dup_list = list(duplicated)
    msg = (
      f"Duplicate column names found: {dup_list}. "
      "Rename or disambiguate columns before writing."
    )
    raise ValueError(msg)

  work = df
  if col_order is not None:
    head = [col for col in col_order if col in work.columns]
    tail = sorted(col for col in work.columns if col not in head)
    if drop_unexpected_cols and tail:
      logger.warning("unexpected_columns_dropped", columns=tail)
      work = work.reindex(columns=head, copy=False)
    else:
      work = work.reindex(columns=head + tail, copy=False)
  else:
    work = work.reindex(columns=sorted(work.columns), copy=False)

  if key_cols is None:
    raise ValueError("key_cols must be provided")
  key_cols_list = list(key_cols)
  sort_cols = _resolve_sort_columns(work, key_cols_list)

  if not sort_cols or sort_chunksize is None or len(work) <= sort_chunksize:
    if sort_cols:
      work.sort_values(by=sort_cols, kind="mergesort", inplace=True)

    for column in work.columns:
      series = work[column]
      if ptypes.is_bool_dtype(series):
        work[column] = _normalise_bool(series)
      else:
        work[column] = _normalise_dates(series)

    with tempfile.NamedTemporaryFile(
      "w", encoding=encoding, newline="\n", delete=False, dir=str(out_path.parent)
    ) as handle:
      tmp_path = Path(handle.name)
      work.to_csv(
        handle,
        index=False,
        float_format="%.6g",
        na_rep="",
        chunksize=chunksize,
        sep=sep,
      )
    os.replace(tmp_path, out_path)
  else:
    with tempfile.TemporaryDirectory() as tmpdir:
      tmp_paths: list[Path] = []
      column_names = list(work.columns)
      for start in range(0, len(work), sort_chunksize):
        subset = work.iloc[start : start + sort_chunksize].sort_values(
          by=sort_cols, kind="mergesort"
        )
        for column in subset.columns:
          series = subset[column]
          if ptypes.is_bool_dtype(series):
            subset[column] = _normalise_bool(series)
          else:
            subset[column] = _normalise_dates(series)
        chunk_path = Path(tmpdir) / f"chunk_{start}.csv"
        subset.to_csv(
          chunk_path,
          index=False,
          float_format="%.6g",
          na_rep="",
          sep=sep,
          encoding=encoding,
        )
        tmp_paths.append(chunk_path)

      key_indices = [column_names.index(col) for col in sort_cols] if sort_cols else []
      key_funcs: list[Callable[[str], object]] = []
      for column in sort_cols:
        dtype = work.dtypes[column]
        if ptypes.is_integer_dtype(dtype):
          key_funcs.append(int)
        elif ptypes.is_float_dtype(dtype):
          key_funcs.append(float)
        else:
          key_funcs.append(lambda value: value)

      with tempfile.NamedTemporaryFile(
        "w",
        encoding=encoding,
        newline="\n",
        delete=False,
        dir=str(out_path.parent),
      ) as handle:
        tmp_path = Path(handle.name)
        writer = csv.writer(handle, delimiter=sep, lineterminator="\n")
        writer.writerow(column_names)

        readers: list[tuple[Iterator[list[str]], TextIO]] = []
        for chunk_file in tmp_paths:
          file_handle = chunk_file.open("r", encoding=encoding, newline="")
          reader = csv.reader(file_handle, delimiter=sep)
          next(reader, None)
          readers.append((reader, file_handle))

        heap: list[tuple[tuple[object, ...], int, list[str]]] = []
        for idx, (reader, _) in enumerate(readers):
          try:
            row = next(reader)
          except StopIteration:
            readers[idx][1].close()
            continue
          key = (
            tuple(func(row[key_indices[pos]]) for pos, func in enumerate(key_funcs))
            if key_funcs
            else tuple()
          )
          heapq.heappush(heap, (key, idx, row))

        while heap:
          _, idx, row = heapq.heappop(heap)
          writer.writerow(row)
          try:
            row = cast(list[str], next(readers[idx][0]))
          except StopIteration:
            readers[idx][1].close()
            continue
          key = (
            tuple(func(row[key_indices[pos]]) for pos, func in enumerate(key_funcs))
            if key_funcs
            else tuple()
          )
          heapq.heappush(heap, (key, idx, row))
    os.replace(tmp_path, out_path)

  write_meta_yaml(
    out_path,
    cfg,
    columns=list(work.columns),
    dtypes={column: work.dtypes[column].name for column in work.columns},
  )
  return out_path


def write_csv_chunks_deterministic(
  chunks: Iterable[pd.DataFrame],
  path: str | Path,
  *,
  col_order: Sequence[str] | None = None,
  key_cols: Sequence[str] | None = None,
  chunksize: int | None = None,
  merge_chunksize: int = 1000,
  sort_chunksize: int | None = None,
  sep: str = ",",
  encoding: str = "utf-8-sig",
  cfg: Config | None = None,
  drop_unexpected_cols: bool = False,
) -> Path:
  """Write DataFrame chunks to ``path`` deterministically."""

  if key_cols is None:
    raise ValueError("key_cols must be provided")
  key_cols_list = list(key_cols)

  with tempfile.TemporaryDirectory() as tmpdir:
    tmp_paths: list[Path] = []
    for idx, df_chunk in enumerate(chunks):
      duplicated = df_chunk.columns[df_chunk.columns.duplicated()]
      if not duplicated.empty:
        dup_list = list(duplicated)
        msg = (
          f"Duplicate column names found in chunk {idx}: {dup_list}. "
          "Rename or disambiguate columns before writing."
        )
        raise ValueError(msg)

      tmp_path = Path(tmpdir) / f"chunk_{idx}.csv"
      write_csv_deterministic(
        df_chunk,
        tmp_path,
        col_order=col_order,
        key_cols=key_cols_list,
        chunksize=chunksize,
        sort_chunksize=sort_chunksize,
        sep=sep,
        encoding=encoding,
        cfg=None,
        drop_unexpected_cols=drop_unexpected_cols,
      )
      meta = tmp_path.with_suffix(tmp_path.suffix + ".meta.yaml")
      if meta.exists():
        meta.unlink()
      tmp_paths.append(tmp_path)

    out_path = Path(path)
    if not tmp_paths:
      return write_csv_deterministic(
        pd.DataFrame(),
        out_path,
        col_order=col_order,
        key_cols=key_cols_list,
        chunksize=chunksize,
        sep=sep,
        encoding=encoding,
        cfg=cfg,
        drop_unexpected_cols=drop_unexpected_cols,
      )

    with ExitStack() as stack:
      readers = []
      for tmp_path in tmp_paths:
        handle = stack.enter_context(tmp_path.open("r", encoding=encoding, newline=""))
        reader = pd.read_csv(
          handle,
          sep=sep,
          chunksize=merge_chunksize,
        )
        stack.callback(reader.close)
        readers.append(reader)
      current: list[pd.DataFrame | None] = []
      for reader in readers:
        try:
          current.append(next(reader))
        except StopIteration:
          current.append(None)

      first = next((chunk for chunk in current if chunk is not None), pd.DataFrame())
      columns = list(first.columns)
      dtypes = {col: first.dtypes[col].name for col in columns}
      resolved_sort_cols = _resolve_sort_columns(first, key_cols_list, emit_warning=False)

      def _fmt(value: Any) -> Any:
        if pd.isna(value):
          return ""
        if isinstance(value, float):
          return f"{value:.6g}"
        if isinstance(value, bool | np.bool_):
          return "true" if bool(value) else "false"
        return value

      heap: list[tuple[tuple[Any, ...], int, int]] = []
      for idx, frame in enumerate(current):
        if frame is not None and not frame.empty:
          row = frame.iloc[0]
          key = (
            tuple(row[col] for col in resolved_sort_cols)
            if resolved_sort_cols
            else tuple()
          )
          heapq.heappush(heap, (key, idx, 0))

      with out_path.open("w", encoding=encoding, newline="") as handle:
        writer = csv.writer(handle, delimiter=sep, lineterminator="\n")
        writer.writerow(columns)
        while heap:
          _, file_idx, row_idx = heapq.heappop(heap)
          chunk = current[file_idx]
          assert chunk is not None
          row = chunk.iloc[row_idx]
          writer.writerow([_fmt(row[column]) for column in columns])
          row_idx += 1
          if row_idx < len(chunk):
            next_row = chunk.iloc[row_idx]
            key = (
              tuple(next_row[col] for col in resolved_sort_cols)
              if resolved_sort_cols
              else tuple()
            )
            heapq.heappush(heap, (key, file_idx, row_idx))
          else:
            try:
              next_chunk = next(readers[file_idx])
            except StopIteration:
              current[file_idx] = None
            else:
              current[file_idx] = next_chunk
              if not next_chunk.empty:
                next_row = next_chunk.iloc[0]
                key = (
                  tuple(next_row[col] for col in resolved_sort_cols)
                  if resolved_sort_cols
                  else tuple()
                )
                heapq.heappush(heap, (key, file_idx, 0))

    write_meta_yaml(out_path, cfg, columns=columns, dtypes=dtypes)
    return out_path


def write_csv(
  df: pd.DataFrame,
  path: str | Path,
  *,
  cfg: Config,
  sep: str | None = None,
  encoding: str | None = None,
  key_cols: Iterable[str] | None = None,
  col_order: Iterable[str] | None = None,
  chunksize: int | None = None,
) -> Path:
  """Write ``df`` to ``path`` as CSV and store metadata."""

  sep = sep or cfg.io.csv_sep
  encoding = encoding or cfg.io.csv_encoding
  key_cols_list = list(key_cols) if key_cols is not None else sorted(df.columns)
  missing_keys = [col for col in key_cols_list if col not in df.columns]
  if missing_keys:
    logger.error(
      "missing_key_columns",
      requested=key_cols_list,
      missing=missing_keys,
    )
    raise ValueError(f"Missing key columns: {missing_keys}")
  col_order_list = list(col_order) if col_order is not None else None
  return write_csv_deterministic(
    df,
    path,
    key_cols=key_cols_list,
    col_order=col_order_list,
    chunksize=chunksize,
    sep=sep,
    encoding=encoding,
    cfg=cfg,
  )


def _read_checkpoint(path: Path) -> int:
  """Return number of rows already processed."""

  try:
    with path.open("r", encoding="utf8") as handle:
      data = json.load(handle)
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
  with tmp.open("w", encoding="utf8") as handle:
    json.dump({"rows": rows}, handle)
  tmp.replace(path)


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
) -> int:
  """Copy ``input_path`` to ``output_path`` using chunked I/O."""

  sep = sep or cfg.csv_sep
  encoding = encoding or cfg.csv_encoding

  processed = 0
  header_written = False
  if checkpoint_path is not None:
    processed = _read_checkpoint(checkpoint_path)
    header_written = Path(output_path).exists() and processed > 0
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
    )
    header_written = True
    rows_written += len(chunk)
    if checkpoint_path is not None:
      _write_checkpoint(checkpoint_path, rows_written)
    if limit is not None and rows_written >= limit:
      break

  return rows_written


def write_meta_yaml(
  path: Path | str,
  cfg: Config | None = None,
  *,
  columns: Sequence[str] | None = None,
  dtypes: Mapping[str, str] | None = None,
) -> Path:
  """Write metadata for ``path`` to ``<path>.meta.yaml``."""

  if dtypes is None and columns is not None:
    dtypes = {column: "string" for column in columns}

  meta = {
    "generated_at": datetime.now().isoformat(),
    "git_sha": git_sha(),
    "command": " ".join(sys.argv),
    "columns": list(columns or []),
    "dtypes": dict(dtypes or {}),
    "config": _serialize_paths(cfg.to_dict()) if cfg is not None else {},
  }
  meta_path = Path(f"{path}.meta.yaml")
  with meta_path.open("w", encoding="utf8") as handle:
    yaml.safe_dump(meta, handle, sort_keys=False)
  return meta_path


def sha256_file(path: Path, *, block_size: int = 64 * 1024) -> str:
  """Return SHA-256 hash of ``path``."""

  if not path.is_file():
    raise FileNotFoundError(f"File not found: {path}")

  digest = hashlib.sha256()
  with path.open("rb") as handle:
    for chunk in iter(lambda: handle.read(block_size), b""):
      digest.update(chunk)
  return digest.hexdigest()


__all__ = [
  "process_csv_chunks",
  "sha256_file",
  "write_csv",
  "write_csv_chunks_deterministic",
  "write_csv_deterministic",
  "write_meta_yaml",
]
