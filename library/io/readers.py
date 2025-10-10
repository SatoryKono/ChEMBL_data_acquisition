"""Reader utilities for CSV-based data pipelines."""

from __future__ import annotations

import csv
import locale
from collections.abc import Hashable, Iterable, Iterator, Mapping, Sequence
from importlib import import_module
from pathlib import Path
from typing import TYPE_CHECKING, Any, NamedTuple

import pandas as pd
from pandas.io.parsers import TextFileReader

from ..common.log import logger
from ..config import IoCfg

if TYPE_CHECKING:  # pragma: no cover - import only for typing
    from .. import validation as _validation_module


_validation_mod: _validation_module | None = None


def _get_validation_module() -> _validation_module:
    """Return the lazily imported :mod:`library.validation` module."""

    global _validation_mod
    if _validation_mod is None:
        from .. import validation as _validation_module_runtime

        _validation_mod = _validation_module_runtime
    return _validation_mod


class _EncodingDecodeError(Exception):
    """Wrap :class:`UnicodeDecodeError` with the attempted encoding."""

    def __init__(self, encoding: str, error: UnicodeDecodeError) -> None:
        super().__init__(str(error))
        self.encoding = encoding
        self.error = error


class _MissingColumnError(ValueError):
    """Raised when the expected identifier column is absent."""

    def __init__(
        self, column: str, path: Path, fieldnames: Sequence[str] | None
    ) -> None:
        message = (
            f"column '{column}' not found in {path}; available columns: {fieldnames}"
        )
        super().__init__(message)
        self.column = column
        self.path = path
        self.fieldnames = fieldnames


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
                raise _MissingColumnError(column, path, reader.fieldnames)
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

    sep_candidates: list[str] = []
    seen_seps: set[str] = set()
    _append_candidate(sep, seen_seps, sep_candidates)
    _append_candidate(
        getattr(cfg, "csv_fallback_separators", None), seen_seps, sep_candidates
    )
    if not sep_candidates:
        sep_candidates.append(cfg.csv_sep)

    def _resolve_candidates() -> Iterator[_CollectedId]:
        last_error: UnicodeDecodeError | None = None
        last_missing: _MissingColumnError | None = None
        for sep_candidate in sep_candidates:
            missing_error: _MissingColumnError | None = None
            for candidate in candidates:
                try:
                    collected = _collect_ids_with_encoding(
                        path_obj,
                        column=column,
                        sep=sep_candidate,
                        encoding=candidate,
                        marker_set=marker_set,
                        keep_na_markers=keep_markers,
                    )
                    try:
                        first_item = next(collected)
                    except StopIteration:
                        return iter(())

                    first_peek_item = first_item
                    remaining_items = collected

                    def _with_peek(
                        peek_item: _CollectedId = first_peek_item,
                        rest: Iterator[_CollectedId] = remaining_items,
                    ) -> Iterator[_CollectedId]:
                        yield peek_item
                        yield from rest

                    if sep_candidate != sep:
                        logger.info(
                            "csv_separator_fallback_used",
                            path=str(path_obj),
                            separator=sep_candidate,
                        )
                    return _with_peek()
                except _MissingColumnError as exc:
                    missing_error = exc
                    if last_missing is None:
                        last_missing = exc
                    break
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
            if missing_error is None:
                # Encoding fallbacks exhausted for this separator; try the next separator.
                continue
        if last_missing is not None:
            raise last_missing
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
    schema: pa.DataFrameSchema | type[pa.DataFrameModel] | None = None,
    chunksize: int | None = None,
) -> pd.DataFrame | TextFileReader[pd.DataFrame]:
    """Load a CSV file into a :class:`pandas.DataFrame` with optional schema validation."""

    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding
    path_obj = Path(path)

    try:
        df = pd.read_csv(
            path_obj,
            sep=sep,
            encoding=encoding,
            dtype=dtype,
            na_values=na_values,
            parse_dates=list(parse_dates) if parse_dates is not None else None,
            chunksize=chunksize,
        )
    except (FileNotFoundError, pd.errors.ParserError, UnicodeError) as exc:
        logger.error(
            "read_fail",
            path=str(path_obj),
            encoding=encoding,
            error=str(exc),
        )
        raise SystemExit(1) from exc
    if chunksize is not None:
        return df
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
        _get_validation_module().validate_columns(df, required_columns)
    return df
