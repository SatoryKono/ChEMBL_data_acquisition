"""Helpers for executing chunked data pipelines.

The utilities centralise the boilerplate involved in running chunked fetchers
with optional concurrency while streaming deterministic CSV output. Pipelines
provide the chunk-level fetch function together with the existing
``run_pipeline`` arguments and the helper coordinates the rest.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from concurrent.futures import (
    FIRST_COMPLETED,
    Future,
    ThreadPoolExecutor,
    as_completed,
    wait,
)
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TypeVar

import pandas as pd

from ...cli.pipeline_definition import normalise_definition
from ...cli_utils import run_pipeline
from ...clients import _chunked

Chunk = TypeVar("Chunk", bound=pd.DataFrame)


@dataclass(slots=True)
class ChunkedFetchConfig:
    """Configuration describing how identifiers should be chunked."""

    ids: Iterable[str]
    chunk_size: int
    workers: int = 1
    chunker: Callable[[Iterable[str], int], Iterable[Sequence[str]]] = _chunked


@dataclass(slots=True)
class CsvWriterConfig:
    """Parameters forwarded to deterministic CSV writers."""

    writer: Callable[..., Path]
    kwargs: Mapping[str, Any]
    ensure_destination: bool = False


RunPipelineKwargs = Mapping[str, Any]


def _ordered_results(
    *,
    chunk_iter: Iterable[Sequence[str]],
    fetch_chunk: Callable[[Sequence[str]], Chunk],
    workers: int,
) -> Iterator[Chunk]:
    """Yield fetched chunks preserving submission order."""

    if workers <= 1:
        for chunk_ids in chunk_iter:
            yield fetch_chunk(list(chunk_ids))
        return

    pending: dict[Future[Chunk], int] = {}
    completed: dict[int, Chunk] = {}
    next_index = 0

    def _cancel_pending() -> None:
        for future in list(pending):
            future.cancel()

    executor = ThreadPoolExecutor(max_workers=workers)
    shutdown_called = False

    def _shutdown_executor(*, wait: bool, cancel_futures: bool) -> None:
        nonlocal shutdown_called
        if shutdown_called:
            return
        shutdown_called = True
        executor.shutdown(wait=wait, cancel_futures=cancel_futures)

    try:
        for index, chunk_ids in enumerate(chunk_iter):
            future = executor.submit(fetch_chunk, list(chunk_ids))
            pending[future] = index
            if len(pending) < workers:
                continue
            done, _ = wait(set(pending), return_when=FIRST_COMPLETED)
            for finished in done:
                chunk_index = pending.pop(finished)
                try:
                    completed[chunk_index] = finished.result()
                except BaseException:
                    _shutdown_executor(wait=False, cancel_futures=True)
                    _cancel_pending()
                    raise
            while next_index in completed:
                yield completed.pop(next_index)
                next_index += 1

        for future in as_completed(list(pending)):
            chunk_index = pending.pop(future)
            try:
                completed[chunk_index] = future.result()
            except BaseException:
                _shutdown_executor(wait=False, cancel_futures=True)
                _cancel_pending()
                raise
            while next_index in completed:
                yield completed.pop(next_index)
                next_index += 1
    finally:
        pending.clear()
        completed.clear()
        _shutdown_executor(wait=True, cancel_futures=False)


def prepare_chunked_pipeline(
    *,
    fetch_config: ChunkedFetchConfig,
    fetch_chunk: Callable[[Sequence[str]], Chunk],
    csv_writer: CsvWriterConfig,
) -> tuple[
    Callable[[], Iterator[Chunk]],
    Callable[[Iterable[Chunk], Path, Sequence[str] | None, Sequence[str]], Path],
]:
    """Construct fetcher and writer callables handling chunking and ordering."""

    def fetcher() -> Iterator[Chunk]:
        chunk_iterable = fetch_config.chunker(
            fetch_config.ids,
            fetch_config.chunk_size,
        )
        yield from _ordered_results(
            chunk_iter=chunk_iterable,
            fetch_chunk=fetch_chunk,
            workers=max(1, fetch_config.workers),
        )

    def writer(
        chunks: Iterable[Chunk],
        destination: Path,
        col_order: Sequence[str] | None,
        key_cols: Sequence[str],
    ) -> Path:
        sort_columns = list(key_cols)
        if not sort_columns and col_order is not None:
            sort_columns = list(col_order)
        order = list(col_order) if col_order is not None else sorted(sort_columns)
        path = csv_writer.writer(
            chunks,
            destination,
            key_cols=sort_columns,
            col_order=order,
            **csv_writer.kwargs,
        )
        path_obj = Path(path)
        if csv_writer.ensure_destination and not path_obj.exists():
            path_obj.parent.mkdir(parents=True, exist_ok=True)
            path_obj.touch()
        return path_obj

    return fetcher, writer


def run_chunked_pipeline(
    *,
    fetch_config: ChunkedFetchConfig,
    fetch_chunk: Callable[[Sequence[str]], Chunk],
    csv_writer: CsvWriterConfig,
    pipeline_kwargs: RunPipelineKwargs,
) -> int:
    """Execute a pipeline using chunked fetching and deterministic CSV output."""

    fetcher, writer = prepare_chunked_pipeline(
        fetch_config=fetch_config,
        fetch_chunk=fetch_chunk,
        csv_writer=csv_writer,
    )

    params = dict(pipeline_kwargs)
    try:
        output_path = params.pop("output_path")
        failure_path = params.pop("failure_path")
    except KeyError as exc:  # pragma: no cover - defensive validation
        missing = exc.args[0]
        raise TypeError(f"run_pipeline missing required argument: {missing}") from exc

    cfg = params.pop("cfg", None)
    logger = params.pop("logger", None)
    definition = params.pop("definition", None)
    emit_standard_outputs = params.pop("emit_standard_outputs", True)
    emit_legacy_artifacts = params.pop("emit_legacy_artifacts", True)
    params.setdefault("writer", writer)

    pipeline_definition = normalise_definition(definition, params)

    return run_pipeline(
        definition=pipeline_definition,
        fetcher=fetcher,
        output_path=output_path,
        failure_path=failure_path,
        cfg=cfg,
        logger=logger,
        emit_standard_outputs=emit_standard_outputs,
        emit_legacy_artifacts=emit_legacy_artifacts,
    )
