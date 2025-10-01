"""Reader utilities for CSV-based data pipelines."""

from __future__ import annotations

import csv
import locale
from collections.abc import Hashable, Iterable, Iterator, Mapping, Sequence
from importlib import import_module
from pathlib import Path
from typing import TYPE_CHECKING, Any, NamedTuple

import pandas as pd

from .. import validation
from ..config import IoCfg
from ..log import logger


class _EncodingDecodeError(Exception):
    """Wrap :class:`UnicodeDecodeError` with the attempted encoding."""

    def __init__(self, encoding: str, error: UnicodeDecodeError) -> None:
        super().__init__(str(error))
        self.encoding = encoding
        self.error = error


class _CollectedId(NamedTuple):
    """Result of parsing a single identifier row."""

    value: str
    is_na_marker: bool
    include: bool


def _collect_ids_with_encoding(
    path: Path,
    *,
    column: str,
    sep: str,
    encoding: str,
    marker_set: set[str],
    keep_na_markers: bool,
) -> Iterator[_CollectedId]:
    """Yield identifiers and NA-marker hits using ``encoding``."""

    try:
        with path.open("r", encoding=encoding, newline="") as fh:
            reader = csv.DictReader(fh, delimiter=sep)
            if reader.fieldnames is None or column not in reader.fieldnames:
                raise ValueError(
                    f"column '{column}' not found in {path}; available columns: {reader.fieldnames}"
                )
            for row in reader:
                value = (row.get(column) or "").strip()
                if not value:
                    continue
                if value in marker_set:
                    include = keep_na_markers
                    yield _CollectedId(value=value, is_na_marker=True, include=include)
                else:
                    yield _CollectedId(value=value, is_na_marker=False, include=True)
    except UnicodeDecodeError as exc:  # pragma: no cover - exercised via fallback tests
        raise _EncodingDecodeError(encoding, exc) from exc


if TYPE_CHECKING:  # pragma: no cover - only for type checking
    import pandera as pa
else:  # pragma: no cover - exercised in tests via monkeypatch
    try:
        pa = import_module("pandera")
    except (ImportError, TypeError):
        pa = None


def read_ids(
    path: str | Path,
    *,
    column: str,
    cfg: IoCfg,
    sep: str | None = None,
    encoding: str | None = None,
    na_markers: Sequence[str] | None = None,
    keep_na_markers: bool | None = None,
) -> Iterator[str]:
    """Yield identifier values from ``column`` in ``path``."""

    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding
    marker_set = set(na_markers or cfg.na_markers or ())
    keep_markers = cfg.keep_na_markers if keep_na_markers is None else keep_na_markers

    def _append_candidate(
        values: Sequence[str] | str | None, seen: set[str], out: list[str]
    ) -> None:
        if values is None:
            return
        if isinstance(values, str):
            iterable: Sequence[str] = (values,)
        else:
            iterable = values
        for value in iterable:
            if not value:
                continue
            key = value.lower()
            if key in seen:
                continue
            out.append(value)
            seen.add(key)

    candidates: list[str] = []
    seen_candidates: set[str] = set()
    _append_candidate((encoding,) if encoding else None, seen_candidates, candidates)
    fallbacks = getattr(cfg, "csv_fallback_encodings", ()) or ()
    _append_candidate(fallbacks, seen_candidates, candidates)

    locale_encoding = locale.getpreferredencoding(False)
    _append_candidate(locale_encoding, seen_candidates, candidates)

    if not candidates:
        candidates.append("utf-8")
        seen_candidates.add("utf-8")

    path_obj = Path(path)

    def _resolve_candidates() -> Iterator[_CollectedId]:
        last_error: UnicodeDecodeError | None = None
        for candidate in candidates:
            try:
                collected = _collect_ids_with_encoding(
                    path_obj,
                    column=column,
                    sep=sep,
                    encoding=candidate,
                    marker_set=marker_set,
                    keep_na_markers=keep_markers,
                )
                try:
                    first_item = next(collected)
                except StopIteration:
                    return iter(())

                def _with_peek() -> Iterator[_CollectedId]:
                    yield first_item
                    yield from collected

                return _with_peek()
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
        collected_iter = _resolve_candidates()
    except FileNotFoundError:
        raise
    except csv.Error as exc:
        raise ValueError(f"malformed CSV in file: {path}: {exc}") from exc

    def _stream() -> Iterator[str]:
        dropped_total = 0
        dropped_unique: list[str] = []
        seen_dropped: set[str] = set()
        try:
            for item in collected_iter:
                if item.is_na_marker:
                    if item.include:
                        yield item.value
                    else:
                        dropped_total += 1
                        if item.value not in seen_dropped:
                            seen_dropped.add(item.value)
                            dropped_unique.append(item.value)
                        continue
                else:
                    yield item.value
        finally:
            if not keep_markers and dropped_total:
                logger.warning(
                    "read_ids_dropped_na_markers",
                    path=str(path_obj),
                    column=column,
                    dropped_total=dropped_total,
                    dropped_ids=dropped_unique,
                )

    return _stream()


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
    schema: "pa.DataFrameSchema" | type["pa.DataFrameModel"] | None = None,
) -> pd.DataFrame:
    """Load a CSV file into a :class:`pandas.DataFrame` with optional schema validation."""

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
        if isinstance(schema, pa.DataFrameSchema):
            df = schema.validate(df)
        else:
            df = schema.to_schema().validate(df)
    elif required_columns is not None:
        validation.validate_columns(df, required_columns)
    return df
