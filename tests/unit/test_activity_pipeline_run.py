"""Unit tests for :mod:`library.pipelines.activity.run`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.pipelines.activity import run as activity_run
from library.pipelines.activity.run import run_activity_pipeline
from library.pipelines.common import (
    ChunkedFetchConfig,
    CsvWriterConfig,
    PipelineRunResult,
)


class _LoggerStub:
    def __init__(self) -> None:
        self.exceptions: list[str] = []

    def exception(
        self, message: str, *args, **kwargs
    ) -> None:  # pragma: no cover - helper
        self.exceptions.append(message % args if args else message)


def _definition_kwargs(tmp_path: Path) -> dict[str, object]:
    return {
        "schema": None,
        "schema_name": "ActivitiesSchema",
        "validators": (),
        "command": "cmd",
        "config_snapshot": {},
        "inputs": {},
        "key_columns": ("activity_id",),
        "table_quality": lambda *_: None,
        "stats_extra": {},
        "stats_callback": None,
        "dictionary_resources": (),
    }


def test_run_activity_pipeline__successful_execution(
    monkeypatch, cfg, tmp_path
) -> None:
    fetch_config = ChunkedFetchConfig(ids=["A1", "A2"], chunk_size=2, workers=1)
    writer_config = CsvWriterConfig(
        writer=lambda *_args, **_kwargs: tmp_path / "out.csv", kwargs={}
    )

    captured: dict[str, object] = {}

    def fake_prepare_chunked_pipeline(*, fetch_config, fetch_chunk, csv_writer):
        captured["fetch_config"] = fetch_config
        captured["csv_writer"] = csv_writer

        def _fetcher():
            return [pd.DataFrame({"activity_id": ["A1"]})]

        def _writer(chunks, destination, col_order, key_cols):
            captured["writer_destination"] = destination
            return Path(destination)

        return _fetcher, _writer

    monkeypatch.setattr(
        activity_run, "prepare_chunked_pipeline", fake_prepare_chunked_pipeline
    )

    def fake_run_pipeline(*, definition, fetcher, output_path, failure_path, **kwargs):
        captured["definition"] = definition
        captured["fetcher"] = fetcher
        captured["output_path"] = Path(output_path)
        return 0

    monkeypatch.setattr(activity_run, "run_pipeline", fake_run_pipeline)

    save_calls: list[tuple[Path, object]] = []

    def fake_save(path: Path, *, cfg) -> None:
        save_calls.append((Path(path), cfg))

    tracker = activity_run.ChunkFailureTracker()
    monkeypatch.setattr(tracker, "save", fake_save, raising=False)

    result = run_activity_pipeline(
        fetch_config=fetch_config,
        metadata_hooks=[lambda frame: frame.assign(flag=True)],
        fetch_chunk=lambda ids: pd.DataFrame({"activity_id": list(ids)}),
        writer_config=writer_config,
        definition_kwargs=_definition_kwargs(tmp_path),
        cfg=cfg,
        logger=_LoggerStub(),
        output_path=tmp_path / "output.csv",
        failure_path=tmp_path / "failure.csv",
        fetch_failure_path=tmp_path / "fetch_failures.csv",
        chunk_failures=tracker,
    )

    assert isinstance(result, PipelineRunResult)
    assert result.exit_code == 0
    assert result.written is True
    assert result.output_path == tmp_path / "output.csv"
    assert captured["definition"].metadata_hooks
    assert save_calls


def test_run_activity_pipeline__failure_result(monkeypatch, cfg, tmp_path) -> None:
    fetch_config = ChunkedFetchConfig(ids=["A1"], chunk_size=1, workers=1)
    writer_config = CsvWriterConfig(
        writer=lambda *_args, **_kwargs: tmp_path / "out.csv", kwargs={}
    )

    monkeypatch.setattr(
        activity_run,
        "prepare_chunked_pipeline",
        lambda **kwargs: (
            lambda: (),
            lambda *a, **k: Path(k.get("destination", tmp_path / "out.csv")),
        ),
    )

    def fake_run_pipeline(**_kwargs):
        return 2

    monkeypatch.setattr(activity_run, "run_pipeline", fake_run_pipeline)

    tracker = activity_run.ChunkFailureTracker()
    save_calls: list[tuple[Path, object]] = []
    monkeypatch.setattr(
        tracker,
        "save",
        lambda path, *, cfg: save_calls.append((Path(path), cfg)),
        raising=False,
    )

    result = run_activity_pipeline(
        fetch_config=fetch_config,
        metadata_hooks=[],
        fetch_chunk=lambda ids: pd.DataFrame({"activity_id": list(ids)}),
        writer_config=writer_config,
        definition_kwargs=_definition_kwargs(tmp_path),
        cfg=cfg,
        logger=_LoggerStub(),
        output_path=tmp_path / "output.csv",
        failure_path=tmp_path / "failure.csv",
        fetch_failure_path=tmp_path / "fetch_failures.csv",
        chunk_failures=tracker,
    )

    assert result.exit_code == 2
    assert result.reason == "pipeline_failed"
    assert result.written is None
    assert save_calls
