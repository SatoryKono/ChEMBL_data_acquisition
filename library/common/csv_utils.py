"""CSV utilities for deterministic output.

This module provides :func:`write_csv_deterministic` which writes a
:pandas.DataFrame to disk in a reproducible manner.  Column order, row
sorting and value serialisation are normalised to guarantee that repeated
runs produce identical files.
"""

from __future__ import annotations

import csv
import hashlib
import heapq
import tempfile
from collections.abc import Callable, Iterable, Sequence
from datetime import date, datetime
from pathlib import Path
from typing import Any, TextIO, cast

import numpy as np
import pandas as pd
from pandas.api import types as ptypes

from ..config import Config
from ..utils.atomic import open_atomic, robust_replace
from .log import logger
from .metadata_writer import write_meta_yaml
from .run_context import get_current


def _metadata_generated_at() -> str | None:
    """Return the active deterministic ``generated_at`` timestamp."""

    context = get_current()
    if context is not None and context.generated_at:
        return context.generated_at
    return None


def _normalise_bool(series: pd.Series[Any]) -> pd.Series[str]:
    """Return ``series`` with booleans converted to ``"true"``/``"false"``.

    ``pandas`` has two boolean dtypes: the classic ``bool`` which cannot
    represent missing values and the nullable ``boolean`` extension type. We
    map both explicit values to lowercase strings and leave missing values as
    ``<NA>`` which :func:`pandas.DataFrame.to_csv` will serialise using
    ``na_rep``.
    """

    return series.map({True: "true", False: "false"}).astype("string")


def _normalise_dates(series: pd.Series[Any]) -> pd.Series[Any]:
    """Return ``series`` with dates formatted as ``YYYY-MM-DD``.

    Both native ``datetime`` columns and object columns containing
    :class:`datetime.date` or :class:`datetime.datetime` values are
    supported. Non-date values are left untouched.
    """

    if ptypes.is_datetime64_any_dtype(series):
        return series.dt.strftime("%Y-%m-%d")

    if (
        ptypes.is_object_dtype(series)
        and series.dropna().map(lambda x: isinstance(x, date | datetime)).all()
    ):
        result: pd.Series[str] = pd.to_datetime(series).dt.strftime("%Y-%m-%d")
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
        msg = "key_cols must contain at least one column"
        raise ValueError(msg)

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


def _make_key_converter(dtype: Any) -> Callable[[Any], tuple[int, Any]]:
    """Return a callable that converts string values to sortable keys."""

    def _is_missing(value: Any) -> bool:
        if pd.isna(value):
            return True
        if isinstance(value, str):
            return value == ""
        return False

    if dtype is not None and ptypes.is_integer_dtype(dtype):
        return lambda value: (1, 0) if _is_missing(value) else (0, int(value))

    if dtype is not None and ptypes.is_float_dtype(dtype):
        return lambda value: (1, 0.0) if _is_missing(value) else (0, float(value))

    return lambda value: (1, "") if _is_missing(value) else (0, value)


class _CSVChunkStream:
    """Lazy reader that streams DataFrame chunks from a CSV file."""

    __slots__ = ("path", "encoding", "sep", "merge_chunksize", "_reader")

    def __init__(
        self,
        path: Path,
        *,
        encoding: str,
        sep: str,
        merge_chunksize: int,
    ) -> None:
        self.path = path
        self.encoding = encoding
        self.sep = sep
        self.merge_chunksize = merge_chunksize
        self._reader: pd.io.parsers.TextFileReader | None = None

    def next(self) -> pd.DataFrame | None:
        """Return the next chunk from the CSV file or ``None`` when exhausted."""

        if self._reader is None:
            self._reader = pd.read_csv(
                self.path,
                sep=self.sep,
                chunksize=self.merge_chunksize,
                dtype="string",
                encoding=self.encoding,
            )
        assert self._reader is not None
        try:
            return next(self._reader)
        except StopIteration:
            self.close()
            return None

    def close(self) -> None:
        """Release the underlying file handle."""

        if self._reader is not None:
            self._reader.close()
            self._reader = None


def _batched_paths(paths: Sequence[Path], size: int) -> Iterable[list[Path]]:
    """Yield ``paths`` in batches of ``size`` preserving order."""

    if size <= 0:
        raise ValueError("merge_file_window must be a positive integer")
    for start in range(0, len(paths), size):
        yield list(paths[start : start + size])


def _merge_sorted_csv_group(
    input_paths: Sequence[Path],
    *,
    destination: Path,
    encoding: str,
    sep: str,
    merge_chunksize: int,
    key_cols: Sequence[str],
    dtype_lookup: dict[str, Any],
    dtype_names_override: dict[str, str] | None,
    fallback_columns: Sequence[str] | None,
) -> tuple[list[str], dict[str, str]]:
    """Merge ``input_paths`` into ``destination`` preserving deterministic order."""

    streams = [
        _CSVChunkStream(
            path,
            encoding=encoding,
            sep=sep,
            merge_chunksize=merge_chunksize,
        )
        for path in input_paths
    ]
    current: list[pd.DataFrame | None] = []
    first_non_empty: pd.DataFrame | None = None
    empty_columns: list[str] | None = None

    try:
        for stream in streams:
            chunk = stream.next()
            current.append(chunk)
            if chunk is None:
                continue
            if not chunk.empty and first_non_empty is None:
                first_non_empty = chunk
            if chunk.empty and empty_columns is None:
                empty_columns = list(chunk.columns)

        if first_non_empty is None:
            if fallback_columns is not None:
                columns = list(fallback_columns)
            elif empty_columns is not None:
                columns = empty_columns
            elif input_paths:
                preview = pd.read_csv(
                    input_paths[0],
                    sep=sep,
                    nrows=0,
                    dtype="string",
                    encoding=encoding,
                )
                columns = list(preview.columns)
            else:
                columns = []
            resolved_sort_cols: list[str] = []
            dtype_names = {
                col: (dtype_names_override or {}).get(col, "string") for col in columns
            }
            with open_atomic(destination, encoding=encoding, newline="") as fh:
                writer = csv.writer(fh, delimiter=sep, lineterminator="\n")
                writer.writerow(columns)
            return columns, dtype_names

        columns = list(first_non_empty.columns)
        resolved_sort_cols = _resolve_sort_columns(
            first_non_empty, key_cols, emit_warning=False
        )
        dtype_names = {
            col: (dtype_names_override or {}).get(col, first_non_empty.dtypes[col].name)
            for col in columns
        }
        key_converters = [
            _make_key_converter(dtype_lookup.get(col)) for col in resolved_sort_cols
        ]

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
                    tuple(
                        converter(row[col])
                        for converter, col in zip(
                            key_converters, resolved_sort_cols, strict=False
                        )
                    )
                    if resolved_sort_cols
                    else tuple()
                )
                heapq.heappush(heap, (key, idx, 0))

        with open_atomic(destination, encoding=encoding, newline="") as fh:
            writer = csv.writer(fh, delimiter=sep, lineterminator="\n")
            writer.writerow(columns)
            while heap:
                _, file_idx, row_idx = heapq.heappop(heap)
                chunk = current[file_idx]
                assert chunk is not None
                row = chunk.iloc[row_idx]
                writer.writerow([_fmt(row[c]) for c in columns])
                row_idx += 1
                if row_idx < len(chunk):
                    next_row = chunk.iloc[row_idx]
                    key = (
                        tuple(
                            converter(next_row[col])
                            for converter, col in zip(
                                key_converters, resolved_sort_cols, strict=False
                            )
                        )
                        if resolved_sort_cols
                        else tuple()
                    )
                    heapq.heappush(heap, (key, file_idx, row_idx))
                else:
                    next_chunk = streams[file_idx].next()
                    current[file_idx] = next_chunk
                    if next_chunk is not None and not next_chunk.empty:
                        next_row = next_chunk.iloc[0]
                        key = (
                            tuple(
                                converter(next_row[col])
                                for converter, col in zip(
                                    key_converters, resolved_sort_cols, strict=False
                                )
                            )
                            if resolved_sort_cols
                            else tuple()
                        )
                        heapq.heappush(heap, (key, file_idx, 0))
        return columns, dtype_names
    finally:
        for stream in streams:
            stream.close()


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
    """Serialise ``df`` to ``path`` as a deterministic CSV file.

    Column names must be unique; duplicate labels will raise a
    :class:`ValueError`.

    Parameters
    ----------
    df:
        DataFrame to be written.
    path:
        Destination file path.
    col_order:
        Optional preferred column ordering. Columns not listed here are
        appended in lexicographical order unless ``drop_unexpected_cols`` is
        ``True``, in which case they are omitted with a warning.
    key_cols:
        Sequence of column names defining row ordering. ``None`` is invalid
        and results in a :class:`ValueError`. Missing columns trigger a
        deterministic fallback that appends remaining columns in their current
        order. An empty sequence is only permitted when ``df`` is empty.
    chunksize:
        Optional number of rows to write per chunk. Passing a value enables
        streaming output via :meth:`pandas.DataFrame.to_csv`, reducing peak
        memory usage for large tables. ``None`` (the default) writes the
        dataframe in a single call.
    sort_chunksize:
        Optional number of rows sorted per chunk. When provided, the
        dataframe is sorted using an external merge-sort algorithm that
        writes intermediate sorted chunks to disk, reducing peak memory
        usage during sorting. ``None`` (the default) performs an in-memory
        sort.
    drop_unexpected_cols:
        When ``True`` and ``col_order`` is provided, columns not present in
        ``col_order`` are dropped from the output with a log warning.

    Notes
    -----
    The input ``df`` is mutated in-place (column order, sorting and normalisation)
    to minimise memory usage. Callers that need to preserve the original data
    should pass ``df.copy()``.

    Returns
    -------
    pathlib.Path
        Path to the written CSV file.

    Raises
    ------
    ValueError
        If ``df`` contains duplicate column names, ``chunksize`` is not a
        positive integer, ``key_cols`` is ``None`` or required sort keys are
        missing.
    """

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Validate optional chunksize to avoid infinite loops or crashes
    if chunksize is not None and chunksize <= 0:
        msg = f"chunksize must be a positive integer, got {chunksize}"
        raise ValueError(msg)
    if sort_chunksize is not None and sort_chunksize <= 0:
        msg = f"sort_chunksize must be a positive integer, got {sort_chunksize}"
        raise ValueError(msg)

    # Ensure column names are unique to avoid ambiguous output
    duplicated = df.columns[df.columns.duplicated()]
    if not duplicated.empty:
        dup_list = list(duplicated)
        msg = (
            f"Duplicate column names found: {dup_list}. "
            "Rename or disambiguate columns before writing."
        )
        raise ValueError(msg)

    # Operations below mutate ``df`` directly to avoid unnecessary copies.
    # Callers should pass ``df.copy()`` if the original frame must remain
    # unchanged.
    work = df

    # Determine column order using ``reindex`` which can avoid materialising
    # a full copy when ``copy=False`` is passed.  This keeps memory overhead
    # proportional to the number of columns rather than the entire DataFrame.
    if col_order is not None:
        head = [c for c in col_order if c in work.columns]
        tail = sorted(c for c in work.columns if c not in head)
        if drop_unexpected_cols and tail:
            logger.warning("unexpected_columns_dropped", columns=tail)
            work = work.reindex(columns=head, copy=False)
        else:
            work = work.reindex(columns=head + tail, copy=False)
    else:
        work = work.reindex(columns=sorted(work.columns), copy=False)

    # Sort rows deterministically using a stable algorithm.
    if key_cols is None:
        raise ValueError("key_cols must be provided")
    key_cols_list = list(key_cols)

    sort_cols = _resolve_sort_columns(work, key_cols_list)

    # In-memory sort when ``sort_chunksize`` is not specified or large enough.
    if not sort_cols or sort_chunksize is None or len(work) <= sort_chunksize:
        if sort_cols:
            work.sort_values(by=sort_cols, kind="mergesort", inplace=True)

        # Normalise bool and date columns without creating intermediary frames.
        for col in work.columns:
            s = work[col]
            if ptypes.is_bool_dtype(s):
                work[col] = _normalise_bool(s)
            else:
                work[col] = _normalise_dates(s)

        # Write via temporary file for atomicity
        with tempfile.NamedTemporaryFile(
            "w", encoding=encoding, newline="\n", delete=False, dir=str(out_path.parent)
        ) as fh:
            tmp_path = Path(fh.name)
            work.to_csv(
                fh,
                index=False,
                float_format="%.6g",
                na_rep="",
                chunksize=chunksize,
                sep=sep,
                lineterminator="\n",
            )
        robust_replace(tmp_path, out_path)
    else:
        # External merge sort writing intermediate chunks to disk.
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_paths: list[Path] = []
            n = len(work)
            column_names = list(work.columns)
            for start in range(0, n, sort_chunksize):
                sub = work.iloc[start : start + sort_chunksize].sort_values(
                    by=sort_cols, kind="mergesort"
                )
                for col in sub.columns:
                    s = sub[col]
                    if ptypes.is_bool_dtype(s):
                        sub[col] = _normalise_bool(s)
                    else:
                        sub[col] = _normalise_dates(s)
                chunk_path = Path(tmpdir) / f"chunk_{start}.csv"
                sub.to_csv(
                    chunk_path,
                    index=False,
                    float_format="%.6g",
                    na_rep="",
                    sep=sep,
                    encoding=encoding,
                    lineterminator="\n",
                )
                tmp_paths.append(chunk_path)

            key_indices = (
                [column_names.index(c) for c in sort_cols] if sort_cols else []
            )
            key_funcs: list[Callable[[str], object]] = []
            for col in sort_cols:
                dtype = work.dtypes[col]
                if ptypes.is_integer_dtype(dtype):
                    key_funcs.append(int)
                elif ptypes.is_float_dtype(dtype):
                    key_funcs.append(float)
                else:
                    key_funcs.append(lambda x: x)

            with tempfile.NamedTemporaryFile(
                "w",
                encoding=encoding,
                newline="\n",
                delete=False,
                dir=str(out_path.parent),
            ) as fh:
                tmp_path = Path(fh.name)
                writer = csv.writer(fh, delimiter=sep, lineterminator="\n")
                writer.writerow(column_names)

                readers: list[tuple[Any, TextIO]] = []
                for p in tmp_paths:
                    f = p.open("r", encoding=encoding, newline="")
                    r = csv.reader(f, delimiter=sep)
                    next(r, None)
                    readers.append((r, f))

                heap: list[tuple[tuple[object, ...], int, list[str]]] = []
                for idx, (reader, _) in enumerate(readers):
                    try:
                        row: list[str] = next(reader)
                    except StopIteration:
                        readers[idx][1].close()
                        continue
                    key = (
                        tuple(
                            func(row[key_indices[i]])
                            for i, func in enumerate(key_funcs)
                        )
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
                        tuple(
                            func(row[key_indices[i]])
                            for i, func in enumerate(key_funcs)
                        )
                        if key_funcs
                        else tuple()
                    )
                    heapq.heappush(heap, (key, idx, row))

            robust_replace(tmp_path, out_path)

    write_meta_yaml(
        out_path,
        config=cfg,
        columns=list(work.columns),
        dtypes={col: work.dtypes[col].name for col in work.columns},
        generated_at=_metadata_generated_at(),
    )
    return out_path


def write_csv_chunks_deterministic(
    chunks: Iterable[pd.DataFrame],
    path: str | Path,
    *,
    col_order: Sequence[str] | None = None,
    key_cols: Sequence[str] | None = None,
    chunksize: int = 1000,
    merge_chunksize: int = 1000,
    sort_chunksize: int | None = None,
    merge_file_window: int = 64,
    sep: str = ",",
    encoding: str = "utf-8-sig",
    cfg: Config | None = None,
    drop_unexpected_cols: bool = False,
) -> Path:
    """Write DataFrame chunks to ``path`` deterministically.

    Column names in each chunk must be unique; duplicates result in a
    :class:`ValueError`.

    Each ``df_chunk`` is immediately normalised and written to a temporary
    file using :func:`write_csv_deterministic` with ``chunksize`` to reduce
    peak memory usage. After all chunks are processed the temporary files are
    reloaded via a :class:`~collections.abc.Generator` and combined with
    :func:`pandas.concat`. The concatenated frame is finally serialised
    deterministically to ``path``. Temporary CSV chunks are re-read with
    ``dtype="string"`` so that textual representations (for example leading
    zeros) are preserved during the merge step.

    Parameters
    ----------
    chunks:
        Iterable yielding :class:`pandas.DataFrame` objects.
    path:
        Destination CSV file.
    col_order:
        Optional preferred column order. Forwarded to
        :func:`write_csv_deterministic`.
    key_cols:
        Columns defining row order. ``None`` is invalid and results in a
        :class:`ValueError`. Missing columns trigger the same deterministic
        fallback used by :func:`write_csv_deterministic`.
    chunksize:
        Number of rows written per chunk.

    merge_chunksize:
        Number of rows loaded from each temporary file during the merge.
        Higher values may improve throughput at the expense of memory.

    merge_file_window:
        Maximum number of temporary chunk files opened simultaneously during
        the merge step. Larger windows increase memory and file descriptor
        usage but reduce the number of intermediate passes.

    sort_chunksize:
        Optional number of rows sorted per chunk during final write. Passed to
        :func:`write_csv_deterministic`.

    sep:
        Field delimiter.
    encoding:
        Output character encoding.
    cfg:
        Optional configuration used for metadata sidecar creation.
    drop_unexpected_cols:
        When ``True`` unexpected columns are dropped.

    Returns
    -------
    pathlib.Path
        Path to the written CSV file.

    Raises
    ------
    ValueError
        If any chunk contains duplicate column names or ``key_cols`` is ``None``.
    """

    if key_cols is None:
        raise ValueError("key_cols must be provided")
    key_cols_list = list(key_cols)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_paths: list[Path] = []
        original_dtypes: dict[str, Any] | None = None
        meta_dtypes: dict[str, str] | None = None
        meta_columns: list[str] | None = None

        for idx, df_chunk in enumerate(chunks):
            # Detect duplicate column names within each chunk
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
            if original_dtypes is None:
                original_dtypes = df_chunk.dtypes.to_dict()
                meta_dtypes = {
                    col: dtype.name if hasattr(dtype, "name") else str(dtype)
                    for col, dtype in original_dtypes.items()
                }
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

        dtype_lookup = original_dtypes or {}
        if merge_file_window <= 0:
            raise ValueError("merge_file_window must be a positive integer")
        pending = list(tmp_paths)
        columns: list[str] = []
        dtype_names: dict[str, str] = {}

        iteration = 0
        while len(pending) > merge_file_window:
            next_round: list[Path] = []
            for batch_index, batch in enumerate(
                _batched_paths(pending, merge_file_window)
            ):
                merged_path = Path(tmpdir) / f"merge_{iteration}_{batch_index}.csv"
                _merge_sorted_csv_group(
                    batch,
                    destination=merged_path,
                    encoding=encoding,
                    sep=sep,
                    merge_chunksize=merge_chunksize,
                    key_cols=key_cols_list,
                    dtype_lookup=dtype_lookup,
                    dtype_names_override=meta_dtypes,
                    fallback_columns=meta_columns,
                )
                next_round.append(merged_path)
                for path_to_remove in batch:
                    if path_to_remove.exists():
                        path_to_remove.unlink()
            pending = next_round
            iteration += 1

        columns, dtype_names = _merge_sorted_csv_group(
            pending,
            destination=out_path,
            encoding=encoding,
            sep=sep,
            merge_chunksize=merge_chunksize,
            key_cols=key_cols_list,
            dtype_lookup=dtype_lookup,
            dtype_names_override=meta_dtypes,
            fallback_columns=meta_columns,
        )

        if meta_columns is None:
            meta_columns = list(columns)

    write_meta_yaml(
        out_path,
        config=cfg,
        columns=meta_columns or columns,
        dtypes=dtype_names,
        generated_at=_metadata_generated_at(),
    )
    return out_path


def sha256_file(path: Path, *, block_size: int = 64 * 1024) -> str:
    """Return SHA-256 hash of ``path``.

    The file is read in chunks to avoid loading large files entirely into
    memory.

    Parameters
    ----------
    path:
        Path to the input file.
    block_size:
        Number of bytes read per iteration. Defaults to ``64 * 1024``.

    Returns
    -------
    str
        Hexadecimal SHA-256 digest of the file contents.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    """

    if not path.is_file():
        raise FileNotFoundError(f"File not found: {path}")

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(block_size), b""):
            digest.update(chunk)
    return digest.hexdigest()
