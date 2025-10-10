"""Unit tests for :mod:`scripts.get_assay_data`."""

from __future__ import annotations

import argparse
import importlib
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import requests
import yaml

from library.cli_utils import PipelineExecutionResult, run_pipeline as cli_run_pipeline
from library.io import StandardOutputArtifacts
from library.config import Config
from library.pipelines.assay.chembl_assay import MAX_ASSAY_CHUNK_SIZE
from library.resources.dictionaries import get_resource
from library.schemas import AssaysSchema
from scripts import get_assay_data


class _MemoryLogger:
    """Capture structured log events emitted by the assay pipeline."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def info(self, event: str, **payload: object) -> None:
        self.events.append(("info", event, dict(payload)))

    def warning(self, event: str, **payload: object) -> None:
        self.events.append(("warning", event, dict(payload)))

    def error(self, event: str, **payload: object) -> None:
        self.events.append(("error", event, dict(payload)))

    def debug(self, event: str, **payload: object) -> None:
        self.events.append(("debug", event, dict(payload)))

    def exception(self, event: str, **payload: object) -> None:
        self.events.append(("exception", event, dict(payload)))


@pytest.fixture()
def logger_stub(monkeypatch: pytest.MonkeyPatch) -> _MemoryLogger:
    logger = _MemoryLogger()
    monkeypatch.setattr(get_assay_data, "logger", logger)
    return logger


@pytest.mark.unit
def test_wrapper_module__exposes_command_public_api() -> None:
    """The compatibility shim must surface the command module attributes."""

    command_module = importlib.import_module("library.cli.commands.get_assay_data")

    assert get_assay_data is command_module
    assert get_assay_data.__all__ == command_module.__all__

    for symbol in [
        "MAX_ASSAY_CHUNK_SIZE",
        "ASSAY_MAX_IDS_PER_REQUEST",
        "_ASSAY_MAX_IDS_PER_REQUEST",
        "ASSAY_OUTPUT_DROP_COLUMNS",
        "_ASSAY_OUTPUT_DROP_COLUMNS",
    ]:
        assert getattr(get_assay_data, symbol) is getattr(command_module, symbol)


@pytest.mark.unit
def test_legacy_assay_max_ids_constant__matches_chunk_size() -> None:
    """Ensure backwards compatible aliases match the public chunk size limit."""

    assert get_assay_data._ASSAY_MAX_IDS_PER_REQUEST == MAX_ASSAY_CHUNK_SIZE
    assert get_assay_data.ASSAY_MAX_IDS_PER_REQUEST == MAX_ASSAY_CHUNK_SIZE


@pytest.fixture()
def minimal_args(tmp_path: Path) -> argparse.Namespace:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("assay_chembl_id\nCHEMBL1\n", encoding="utf-8")
    final_out = tmp_path / "output.csv"
    return argparse.Namespace(
        input_csv=input_csv,
        final_out=final_out,
        skip_existing=False,
        force=False,
        offset=0,
        emit_legacy_artifacts=False,
    )


def test_run_chembl__invalid_limit_logs_error(
    cfg: Config, minimal_args: argparse.Namespace, logger_stub: _MemoryLogger
) -> None:
    cfg.assay.limit = -1

    exit_code = get_assay_data.run_chembl(cfg, minimal_args)

    assert exit_code == 1
    assert (
        "error",
        "invalid_limit",
        {"section": "assay.limit", "limit": -1},
    ) in logger_stub.events


def test_run_chembl__read_ids_failure(
    cfg: Config,
    minimal_args: argparse.Namespace,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_read_ids(path: Path, *, column: str, cfg: Any) -> Any:
        raise FileNotFoundError(path)

    monkeypatch.setattr(get_assay_data.io, "read_ids", fake_read_ids)

    exit_code = get_assay_data.run_chembl(cfg, minimal_args)

    assert exit_code == 1
    assert any(event == "read_fail" for _, event, _ in logger_stub.events)


@pytest.mark.parametrize("offset", [0, 2])
def test_run_chembl__successful_execution(
    cfg: Config,
    minimal_args: argparse.Namespace,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
    offset: int,
) -> None:
    minimal_args.offset = offset
    cfg.assay.limit = 2

    outer_cfg = cfg

    def fake_read_ids(path: Path, *, column: str, cfg: Any) -> Any:
        assert column == outer_cfg.assay.column
        assert cfg is outer_cfg.io
        return iter(["CHEMBL1", "CHEMBL2", "CHEMBL3"])

    class FakeClient:
        def __enter__(self) -> FakeClient:
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

        def close(self) -> None:
            return None

    class FakeTracker:
        def __init__(self) -> None:
            self.stats: dict[str, object] = {"failures": 0}

        def add_failure(self, *_: object, **__: object) -> None:
            return None

        def save(self, path: Path, *, cfg: Config) -> None:
            self.saved_path = path

    tracker = FakeTracker()

    def fake_prepare_chunked_pipeline(**kwargs: object):
        def fetcher() -> Iterable[pd.DataFrame]:
            yield pd.DataFrame({"assay_chembl_id": ["CHEMBL1"]})

        def writer(**_: object) -> Path:
            return minimal_args.final_out

        return fetcher, writer

    def fake_run_pipeline(
        *,
        definition,
        fetcher,
        output_path,
        failure_path,
        **kwargs: object,
    ) -> PipelineExecutionResult:
        del fetcher, output_path, failure_path, kwargs
        if definition.stats_callback is not None:
            definition.stats_callback(
                {"rows_total": 2, "rows_kept": 2, "rows_dropped": 0}
            )
        pd.DataFrame({"assay_chembl_id": ["CHEMBL1"]}).to_csv(
            minimal_args.final_out,
            index=False,
        )
        return PipelineExecutionResult(
            exit_code=0,
            dataset_path=minimal_args.final_out,
        )

    def fake_save_standard_outputs(
        dataset: pd.DataFrame,
        quality: pd.DataFrame,
        correlation: pd.DataFrame,
        **_: object,
    ) -> StandardOutputArtifacts:
        del dataset, quality, correlation
        return StandardOutputArtifacts(
            dataset=minimal_args.final_out,
            quality_report=minimal_args.final_out.with_name("quality.csv"),
            correlation_report=minimal_args.final_out.with_name(
                "correlation.csv"
            ),
        )

    monkeypatch.setattr(get_assay_data.io, "read_ids", fake_read_ids)
    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        lambda *args, **kwargs: FakeClient(),
    )
    monkeypatch.setattr(get_assay_data, "ChunkFailureTracker", lambda: tracker)
    monkeypatch.setattr(
        get_assay_data.cl,
        "get_assays",
        lambda *args, **kwargs: pd.DataFrame({"assay_chembl_id": ["CHEMBL1"]}),
    )
    monkeypatch.setattr(
        get_assay_data, "prepare_chunked_pipeline", fake_prepare_chunked_pipeline
    )
    monkeypatch.setattr(get_assay_data, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(
        get_assay_data.io,
        "save_standard_outputs",
        fake_save_standard_outputs,
    )

    exit_code = get_assay_data.run_chembl(cfg, minimal_args)

    assert exit_code == 0
    assert any(
        event == "assay_standard_outputs" for _, event, _ in logger_stub.events
    )
    assert any(event == "postprocess_skipped" for _, event, _ in logger_stub.events)
    if offset:
        assert any(event == "process_offset" for _, event, _ in logger_stub.events)
    assert any(event == "process_limit" for _, event, _ in logger_stub.events)


@pytest.mark.unit
def test_run_chembl__splits_chunk_on_timeout(
    cfg: Config,
    minimal_args: argparse.Namespace,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg.assay.limit = None
    cfg.assay.batch_size = 4
    cfg.retry.max_attempts = 2
    cfg.retry.backoff_factor = 0

    def fake_read_ids(path: Path, *, column: str, cfg: Any) -> Iterable[str]:
        return iter(["CHEMBL100", "CHEMBL200"])

    class FakeClient:
        def __enter__(self) -> FakeClient:
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

        def close(self) -> None:
            return None

    class FakeTracker:
        def __init__(self) -> None:
            self._failures: list[tuple[list[str], str]] = []

        def add_failure(self, chunk_ids: Iterable[str], error: str) -> None:
            self._failures.append((list(chunk_ids), error))

        def save(self, path: Path, *, cfg: Config) -> None:
            self.saved_path = path

        def stats(self) -> dict[str, object]:
            return {"failures": len(self._failures)}

    tracker = FakeTracker()
    call_history: list[list[str]] = []

    def fake_get_assays(
        identifiers: Sequence[str],
        *,
        cfg: Any,
        client: Any,
        chunk_size: int,
        timeout: float,
    ) -> pd.DataFrame:
        call_history.append(list(identifiers))
        if len(identifiers) > 1:
            raise requests.ReadTimeout("timeout while fetching chunk")
        return pd.DataFrame({"assay_chembl_id": list(identifiers)})

    def fake_prepare_chunked_pipeline(**kwargs: object):
        fetch_chunk = kwargs["fetch_chunk"]

        def fetcher() -> Iterable[pd.DataFrame]:
            yield fetch_chunk(["CHEMBL100", "CHEMBL200"])

        def writer(
            frames: Iterable[pd.DataFrame],
            destination: Path,
            col_order: Sequence[str] | None,
            key_cols: Sequence[str],
        ) -> Path:
            del destination, col_order, key_cols
            list(frames)
            return minimal_args.final_out

        return fetcher, writer

    def fake_run_pipeline(
        *,
        definition,
        fetcher: Callable[[], Iterable[pd.DataFrame]],
        output_path: Path,
        failure_path: Path,
        **kwargs: object,
    ) -> PipelineExecutionResult:
        del output_path, failure_path, kwargs
        list(fetcher())
        if definition.stats_callback is not None:
            definition.stats_callback({})
        pd.DataFrame({"assay_chembl_id": ["CHEMBL100", "CHEMBL200"]}).to_csv(
            minimal_args.final_out,
            index=False,
        )
        return PipelineExecutionResult(
            exit_code=0,
            dataset_path=minimal_args.final_out,
        )

    def fake_save_standard_outputs(
        dataset: pd.DataFrame,
        quality: pd.DataFrame,
        correlation: pd.DataFrame,
        **_: object,
    ) -> StandardOutputArtifacts:
        return StandardOutputArtifacts(
            dataset=minimal_args.final_out,
            quality_report=minimal_args.final_out.with_name("quality.csv"),
            correlation_report=minimal_args.final_out.with_name(
                "correlation.csv"
            ),
        )

    monkeypatch.setattr(get_assay_data.io, "read_ids", fake_read_ids)
    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        lambda *_, **__: FakeClient(),
    )
    monkeypatch.setattr(get_assay_data, "ChunkFailureTracker", lambda: tracker)
    monkeypatch.setattr(get_assay_data.cl, "get_assays", fake_get_assays)
    monkeypatch.setattr(
        get_assay_data, "prepare_chunked_pipeline", fake_prepare_chunked_pipeline
    )
    monkeypatch.setattr(get_assay_data, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(get_assay_data, "sleep", lambda *_: None)
    monkeypatch.setattr(
        get_assay_data.io,
        "save_standard_outputs",
        fake_save_standard_outputs,
    )

    exit_code = get_assay_data.run_chembl(cfg, minimal_args)

    assert exit_code == 0
    assert call_history[0] == ["CHEMBL100", "CHEMBL200"]
    assert call_history.count(["CHEMBL100", "CHEMBL200"]) == 2
    assert ["CHEMBL100"] in call_history
    assert ["CHEMBL200"] in call_history
    assert tracker._failures == []
    assert any(event == "assay_fetch_split" for _, event, _ in logger_stub.events)


def test_run_chembl__standard_outputs_created_without_legacy(
    cfg: Config,
    minimal_args: argparse.Namespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg.assay.limit = 1

    df = pd.DataFrame({"assay_chembl_id": ["CHEMBL1"]})
    df.to_csv(minimal_args.final_out, index=False, encoding="utf-8")

    save_calls: list[tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]] = []

    def fake_save_standard_outputs(dataset, quality, correlation, **kwargs: object):
        save_calls.append((dataset, quality, correlation))
        return StandardOutputArtifacts(
            dataset=minimal_args.final_out,
            quality_report=minimal_args.final_out.with_name("quality.csv"),
            correlation_report=minimal_args.final_out.with_name("correlation.csv"),
        )

    class FakeTracker:
        def __init__(self) -> None:
            self.saved = False

        def add_failure(self, *_: object, **__: object) -> None:
            return None

        def save(self, *_: object, **__: object) -> None:
            self.saved = True

        def stats(self) -> dict[str, object]:
            return {"failures": 0}

    tracker = FakeTracker()

    def fake_run_pipeline(**kwargs: object) -> PipelineExecutionResult:
        if (definition := kwargs.get("definition")) and definition.stats_callback:
            definition.stats_callback({"rows_total": 1, "rows_kept": 1, "rows_dropped": 0})
        return PipelineExecutionResult(
            exit_code=0,
            dataset_path=minimal_args.final_out,
        )

    monkeypatch.setattr(get_assay_data, "ChunkFailureTracker", lambda: tracker)
    monkeypatch.setattr(get_assay_data, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(
        get_assay_data.io,
        "save_standard_outputs",
        fake_save_standard_outputs,
    )
    monkeypatch.setattr(
        get_assay_data.io,
        "read_ids",
        lambda *_args, **_kwargs: iter(["CHEMBL1"]),
    )
    monkeypatch.setattr(get_assay_data.cl, "get_assays", lambda *_, **__: df)
    class _DummyClient:
        def __enter__(self) -> "_DummyClient":
            return self

        def __exit__(self, *_: object) -> None:
            return None

        def close(self) -> None:
            return None

    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        lambda *_, **__: _DummyClient(),
    )
    def fake_prepare_chunked_pipeline(**kwargs: object):
        del kwargs

        def _fetcher() -> Iterable[pd.DataFrame]:
            yield df

        def _writer(*_args: object, **_kwargs: object) -> Path:
            return minimal_args.final_out

        return _fetcher, _writer

    monkeypatch.setattr(
        get_assay_data,
        "prepare_chunked_pipeline",
        fake_prepare_chunked_pipeline,
    )

    exit_code = get_assay_data.run_chembl(cfg, minimal_args)

    assert exit_code == 0
    assert save_calls, "expected save_standard_outputs to be invoked"
    assert not tracker.saved
    fetch_failure = minimal_args.final_out.with_name(
        f"{minimal_args.final_out.stem}_fetch_failures.csv"
    )
    assert not fetch_failure.exists()


@pytest.mark.unit
def test_run_chembl__quality_failure_respects_fatal_flag(
    cfg: Config,
    minimal_args: argparse.Namespace,
    monkeypatch: pytest.MonkeyPatch,
    logger_stub: _MemoryLogger,
) -> None:
    cfg.assay.limit = 1
    cfg.system.doc_quality.fatal_on_error = True

    df = pd.DataFrame({"assay_chembl_id": ["CHEMBL1"]})
    df.to_csv(minimal_args.final_out, index=False, encoding="utf-8")

    monkeypatch.setattr(
        get_assay_data.io,
        "read_ids",
        lambda *_args, **_kwargs: iter(["CHEMBL1"]),
    )

    def fake_prepare_chunked_pipeline(**_kwargs: object):
        def fetcher() -> Iterable[pd.DataFrame]:
            yield df

        def writer(**_writer_kwargs: object) -> Path:
            return minimal_args.final_out

        return fetcher, writer

    monkeypatch.setattr(
        get_assay_data, "prepare_chunked_pipeline", fake_prepare_chunked_pipeline
    )

    def fake_run_pipeline(**_kwargs: object) -> PipelineExecutionResult:
        df.to_csv(minimal_args.final_out, index=False, encoding="utf-8")
        return PipelineExecutionResult(exit_code=0, dataset_path=minimal_args.final_out)

    monkeypatch.setattr(get_assay_data, "run_pipeline", fake_run_pipeline)

    def fail_save_standard_outputs(*_args: object, **_kwargs: object) -> None:
        raise AssertionError(
            "standard outputs should not be persisted after fatal QA failure"
        )

    monkeypatch.setattr(
        get_assay_data.io, "save_standard_outputs", fail_save_standard_outputs
    )

    def failing_quality_hook(*_args: object, **_kwargs: object):
        def _hook(_frame: pd.DataFrame) -> None:
            raise RuntimeError("quality hook failed")

        return _hook

    monkeypatch.setattr(
        get_assay_data, "build_table_quality_hook", failing_quality_hook
    )

    exit_code = get_assay_data.run_chembl(cfg, minimal_args)

    assert exit_code == 1
    assert any(
        event == "assay_quality_generation_failed" and payload.get("fatal")
        for _, event, payload in logger_stub.events
    )


def test_run__skip_existing_returns_zero(
    cfg: Config,
    minimal_args: argparse.Namespace,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    minimal_args.skip_existing = True
    minimal_args.force = False
    minimal_args.final_out.write_text("existing", encoding="utf-8")

    called = False

    def fake_run(cfg: Config, args: argparse.Namespace) -> int:
        nonlocal called
        called = True
        return 0

    monkeypatch.setattr(get_assay_data, "run_chembl", fake_run)

    exit_code = get_assay_data.run(cfg, minimal_args)

    assert exit_code == 0
    assert not called
    assert (
        "info",
        "pipeline_skip_existing",
        {"output": str(minimal_args.final_out)},
    ) in logger_stub.events


def test_run__force_overrides_skip(
    cfg: Config,
    minimal_args: argparse.Namespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    minimal_args.skip_existing = True
    minimal_args.force = True
    minimal_args.final_out.write_text("existing", encoding="utf-8")

    calls: list[str] = []

    def fake_run(cfg: Config, args: argparse.Namespace) -> int:
        calls.append("run")
        return 3

    monkeypatch.setattr(get_assay_data, "run_chembl", fake_run)

    exit_code = get_assay_data.run(cfg, minimal_args)

    assert exit_code == 3
    assert calls == ["run"]


def test_run__propagates_exit_code(
    cfg: Config, minimal_args: argparse.Namespace, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(get_assay_data, "run_chembl", lambda *_: 7)

    exit_code = get_assay_data.run(cfg, minimal_args)

    assert exit_code == 7


@pytest.mark.unit
def test_run_pipeline__adds_missing_assay_optional_columns(tmp_path: Path) -> None:
    frame = pd.DataFrame({"assay_chembl_id": ["CHEMBL1"]})

    def fetcher() -> Iterable[pd.DataFrame]:
        yield frame

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: Iterable[str] | None,
        key_cols: Iterable[str] | None,
        **_: object,
    ) -> Path:
        frames = [chunk.copy() for chunk in chunks]
        if frames:
            result = pd.concat(frames, ignore_index=True)
        else:
            result = pd.DataFrame(columns=list(col_order or []))
        if col_order:
            result = result.reindex(columns=list(col_order), fill_value=pd.NA)
        destination.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(destination, index=False)
        return destination

    logger = _MemoryLogger()
    output_path = tmp_path / "assays.csv"
    failure_path = tmp_path / "assays_failures.csv"

    exit_code = cli_run_pipeline(
        fetcher=fetcher,
        schema=AssaysSchema,
        schema_name="AssaysSchema",
        validators=[],
        metadata_hooks=[],
        writer=writer,
        output_path=output_path,
        failure_path=failure_path,
        command="pytest",
        config_snapshot={},
        inputs={},
        key_columns=["assay_chembl_id"],
        table_quality=lambda _: None,
        cfg=None,
        stats_extra=None,
        logger=logger,
        dictionary_resources=("dictionary_root",),
        emit_legacy_artifacts=False,
    )

    assert exit_code == 0
    assert output_path.exists()
    result = pd.read_csv(output_path)
    assert "assay_group" in result.columns
    assert "assay_strain" in result.columns
    assert result["assay_group"].isna().all()
    assert result["assay_strain"].isna().all()

    meta_path = output_path.with_name(output_path.name + ".meta.yaml")
    assert not meta_path.exists()


def test_build_parser__defaults() -> None:
    parser, _ = get_assay_data.build_parser()
    args = parser.parse_args([])

    assert args.input_csv == Path(get_assay_data.DEFAULT_INPUT_NAME)
    assert args.batch_size == parser.get_default("batch_size")
    assert callable(args.func)
