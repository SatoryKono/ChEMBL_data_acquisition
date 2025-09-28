"""CSV reading utilities for data acquisition pipelines."""

from __future__ import annotations

import csv
import locale
from collections.abc import Hashable, Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from .. import validation
from ..config import IoCfg
from ..utils.logging import logger

if TYPE_CHECKING:  # pragma: no cover - only required for typing
  import pandera as pa
else:  # pragma: no cover - exercised indirectly via runtime imports
  try:
    import pandera as pa  # type: ignore[import-not-found]
  except (ImportError, TypeError):
    pa = None  # type: ignore[assignment]


class _EncodingDecodeError(Exception):
  """Wrap :class:`UnicodeDecodeError` with the attempted encoding."""

  def __init__(self, encoding: str, error: UnicodeDecodeError) -> None:
    super().__init__(str(error))
    self.encoding = encoding
    self.error = error


def _stream_ids_with_encoding(
  path: Path,
  *,
  column: str,
  sep: str,
  encoding: str,
  marker_set: set[str],
) -> Iterator[str]:
  """Yield identifier values using a specific ``encoding``."""

  try:
    with path.open("r", encoding=encoding, newline="") as handle:
      reader = csv.DictReader(handle, delimiter=sep)
      if reader.fieldnames is None or column not in reader.fieldnames:
        available = reader.fieldnames or []
        raise ValueError(
          f"column '{column}' not found in {path}; available columns: {available}"
        )
      for row in reader:
        value = (row.get(column) or "").strip()
        if value and value not in marker_set:
          yield value
  except UnicodeDecodeError as exc:  # pragma: no cover - triggered via fallbacks
    raise _EncodingDecodeError(encoding, exc) from exc


def read_ids(
  path: str | Path,
  *,
  column: str,
  cfg: IoCfg,
  sep: str | None = None,
  encoding: str | None = None,
  na_markers: Sequence[str] | None = None,
) -> Iterator[str]:
  """Yield identifier values from ``column`` in ``path``."""

  sep = sep or cfg.csv_sep
  encoding = encoding or cfg.csv_encoding
  marker_set = set(na_markers or cfg.na_markers or ())

  def append_candidate(
    values: Sequence[str] | str | None,
    seen: set[str],
    out: list[str],
  ) -> None:
    if values is None:
      return
    iterable: Sequence[str]
    if isinstance(values, str):
      iterable = (values,)
    else:
      iterable = values
    for candidate in iterable:
      if not candidate:
        continue
      key = candidate.lower()
      if key in seen:
        continue
      out.append(candidate)
      seen.add(key)

  candidates: list[str] = []
  seen_candidates: set[str] = set()
  append_candidate((encoding,) if encoding else None, seen_candidates, candidates)
  fallbacks = getattr(cfg, "csv_fallback_encodings", ()) or ()
  append_candidate(fallbacks, seen_candidates, candidates)

  locale_encoding = locale.getpreferredencoding(False)
  append_candidate(locale_encoding, seen_candidates, candidates)

  if not candidates:
    candidates.append("utf-8")
    seen_candidates.add("utf-8")

  path_obj = Path(path)

  def iter_candidates() -> Iterator[str]:
    last_error: UnicodeDecodeError | None = None
    for candidate in candidates:
      try:
        yield from _stream_ids_with_encoding(
          path_obj,
          column=column,
          sep=sep,
          encoding=candidate,
          marker_set=marker_set,
        )
        return
      except _EncodingDecodeError as exc:
        last_error = exc.error
        logger.warning(
          "csv_decode_failed",
          path=str(path_obj),
          encoding=exc.encoding,
          error=str(exc.error),
        )
        continue
      except LookupError as exc:
        logger.warning(
          "csv_encoding_lookup_failed",
          path=str(path_obj),
          encoding=candidate,
          error=str(exc),
        )
        continue
    attempted = ", ".join(candidates)
    message = (
      f"failed to decode CSV {path_obj} with encodings: {attempted}. "
      "Update 'io.csv_encoding' or 'io.csv_fallback_encodings' in the configuration."
    )
    raise ValueError(message) from last_error

  try:
    yield from iter_candidates()
  except FileNotFoundError:
    raise
  except csv.Error as exc:
    raise ValueError(f"malformed CSV in file: {path}: {exc}") from exc


def read_csv(
  path: str | Path,
  *,
  cfg: IoCfg,
  sep: str | None = None,
  encoding: str | None = None,
  required_columns: Iterable[str] | None = None,
  dtype: Mapping[Hashable, Any] | type | None = None,
  na_values: Sequence[str] | str | None = None,
  parse_dates: Sequence[str] | None = None,
  schema: "pa.DataFrameSchema | type[pa.DataFrameModel] | None" = None,
) -> pd.DataFrame:
  """Load a CSV file into a :class:`pandas.DataFrame` with optional validation."""

  sep = sep or cfg.csv_sep
  encoding = encoding or cfg.csv_encoding
  df = pd.read_csv(
    path,
    sep=sep,
    encoding=encoding,
    dtype=dtype,
    na_values=na_values,
    parse_dates=list(parse_dates) if parse_dates is not None else None,
  )
  if schema is not None:
    if pa is None:
      raise RuntimeError(
        "pandera is required for schema validation; install pandera to use the 'schema' argument"
      )
    if isinstance(schema, pa.DataFrameSchema):  # type: ignore[arg-type]
      df = schema.validate(df)
    else:
      df = schema.to_schema().validate(df)
  elif required_columns is not None:
    validation.validate_columns(df, required_columns)
  return df


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
  """Yield chunked dataframes from ``path`` honouring ``cfg`` defaults."""

  if chunk_size <= 0:
    raise ValueError(f"chunk_size must be positive, got {chunk_size}")

  sep = sep or cfg.csv_sep
  encoding = encoding or cfg.csv_encoding
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
        if remaining <= 0:
          break
        if len(chunk) > remaining:
          yield chunk.iloc[:remaining]
          emitted += remaining
          break
      emitted += len(chunk)
      yield chunk
  finally:
    reader.close()


__all__ = [
  "read_csv",
  "read_csv_chunks",
  "read_ids",
]
