"""Tests for :mod:`library.pipelines.target.pipeline`."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import Any

import pandas as pd
import pytest

from library.pipelines.target.pipeline import PipelineResult, run_pipeline


def _iterator() -> Iterator[Iterable[str]]:
    yield ["CHEMBL1", "CHEMBL2"]


def test_run_pipeline_passes_batch_size() -> None:
    """Ensure ``batch_size`` is forwarded to fetchers that accept it."""

    received: dict[str, Any] = {}

    def chembl_fetcher(
        chunks: Iterator[Iterable[str]],
        cfg: object,
        *,
        batch_size: int,
    ) -> pd.DataFrame:
        received["batch_size"] = batch_size
        ids = [item for chunk in chunks for item in chunk]
        return pd.DataFrame({"target_chembl_id": ids})

    result = run_pipeline(
        _iterator,
        chembl_cfg=object(),
        chembl_fetcher=chembl_fetcher,
        batch_size=25,
    )

    assert isinstance(result, PipelineResult)
    assert received["batch_size"] == 25
    assert result.chembl.loc[0, "target_chembl_id"] == "CHEMBL1"


def test_run_pipeline_ignores_optional_keywords_for_wrappers() -> None:
    """Wrappers without ``batch_size`` support are retried without the keyword."""

    calls = {"count": 0}

    def base_fetcher(chunks: Iterator[Iterable[str]], cfg: object) -> pd.DataFrame:
        calls["count"] += 1
        ids = [item for chunk in chunks for item in chunk]
        return pd.DataFrame({"target_chembl_id": ids})

    def cached_fetcher(chunks: Iterator[Iterable[str]], cfg: object) -> pd.DataFrame:
        # Simulate a closure that intentionally omits ``batch_size``.
        return base_fetcher(chunks, cfg)

    result = run_pipeline(
        _iterator,
        chembl_cfg=object(),
        chembl_fetcher=cached_fetcher,
        batch_size=15,
    )

    assert calls["count"] == 1
    assert list(result.chembl["target_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]


def test_optional_stages_receive_dataframe() -> None:
    """Optional stages operate on the ChEMBL output."""

    def chembl_fetcher(
        chunks: Iterator[Iterable[str]],
        cfg: object,
        *,
        batch_size: int,
    ) -> pd.DataFrame:
        ids = [item for chunk in chunks for item in chunk]
        return pd.DataFrame({"target_chembl_id": ids})

    def uniprot_fetcher(
        df: pd.DataFrame, cfg: object | None = None, **_: Any
    ) -> pd.DataFrame:
        assert list(df["target_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
        return pd.DataFrame({"uniprot_id": ["P1", "P2"]})

    result = run_pipeline(
        _iterator,
        chembl_cfg=object(),
        chembl_fetcher=chembl_fetcher,
        uniprot_fetcher=uniprot_fetcher,
        uniprot_cfg=object(),
    )

    assert result.uniprot is not None
    assert list(result.uniprot["uniprot_id"]) == ["P1", "P2"]

    # Mapping interface exposes the stored frames.
    assert result["chembl"].equals(result.chembl)
    assert result.as_dict()["uniprot"].equals(result.uniprot)


def test_run_pipeline_streams_frames_with_metadata(monkeypatch: pytest.MonkeyPatch) -> None:
    """Streaming results preserve chunk order and inject metadata lazily."""

    call_sequence: list[int] = []

    def marker(frame: pd.DataFrame) -> pd.DataFrame:
        call_index = len(call_sequence) + 1
        call_sequence.append(call_index)
        # Provide deterministic metadata for assertions.
        return frame.assign(pipeline_marker=call_index)

    monkeypatch.setattr(
        "library.pipelines.target.pipeline.add_pipeline_metadata", marker
    )

    def chunk_iterator() -> Iterator[Iterable[str]]:
        yield ["CHEMBL10", "CHEMBL11"]
        yield ["CHEMBL12"]

    def chembl_fetcher(
        chunks: Iterator[Iterable[str]], cfg: object, *, batch_size: int
    ) -> Iterator[pd.DataFrame]:
        assert batch_size == 100  # default propagated value

        def stream() -> Iterator[pd.DataFrame]:
            for chunk in chunks:
                yield pd.DataFrame({"target_chembl_id": list(chunk)})

        return stream()

    result = run_pipeline(
        chunk_iterator,
        chembl_cfg=object(),
        chembl_fetcher=chembl_fetcher,
    )

    assert isinstance(result.chembl, Iterator)

    streamed_frames = list(result.chembl)
    assert len(streamed_frames) == 2
    assert [
        list(frame["target_chembl_id"])
        for frame in streamed_frames
    ] == [["CHEMBL10", "CHEMBL11"], ["CHEMBL12"]]

    for idx, frame in enumerate(streamed_frames, start=1):
        assert set(frame["pipeline_marker"]) == {idx}

    # Generator is exhausted after a single pass.
    assert list(result.chembl) == []
    assert call_sequence == [1, 2]


def test_run_pipeline_streams_empty_frame(monkeypatch: pytest.MonkeyPatch) -> None:
    """When no chunks are emitted a single empty frame is yielded with metadata."""

    call_sequence: list[int] = []

    def marker(frame: pd.DataFrame) -> pd.DataFrame:
        call_index = len(call_sequence) + 1
        call_sequence.append(call_index)
        return frame.assign(pipeline_marker=call_index)

    monkeypatch.setattr(
        "library.pipelines.target.pipeline.add_pipeline_metadata", marker
    )

    def chunk_iterator() -> Iterator[Iterable[str]]:
        return iter(())

    def chembl_fetcher(
        chunks: Iterator[Iterable[str]], cfg: object, *, batch_size: int
    ) -> Iterator[pd.DataFrame]:

        def stream() -> Iterator[pd.DataFrame]:
            for chunk in chunks:
                yield pd.DataFrame({"target_chembl_id": list(chunk)})

        return stream()

    result = run_pipeline(
        chunk_iterator,
        chembl_cfg=object(),
        chembl_fetcher=chembl_fetcher,
    )

    assert isinstance(result.chembl, Iterator)
    streamed_frames = list(result.chembl)
    assert len(streamed_frames) == 1

    empty_frame = streamed_frames[0]
    assert empty_frame.empty
    assert list(empty_frame.columns) == ["pipeline_marker"]
    assert list(result.chembl) == []
    assert call_sequence == [1]
