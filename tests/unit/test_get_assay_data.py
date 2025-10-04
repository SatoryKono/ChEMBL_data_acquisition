"""Unit tests for :mod:`scripts.get_assay_data`."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import pytest

from library.config import Config
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


@pytest.fixture()
def logger_stub(monkeypatch: pytest.MonkeyPatch) -> _MemoryLogger:
    logger = _MemoryLogger()
    monkeypatch.setattr(get_assay_data, "logger", logger)
    return logger


@pytest.fixture()
def minimal_args(tmp_path: Path) -> argparse.Namespace:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("assay_chembl_id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"
    return argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        skip_existing=False,
        force=False,
        offset=0,
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


def test_run_chembl__read_ids_failure(cfg: Config, minimal_args: argparse.Namespace, logger_stub: _MemoryLogger, monkeypatch: pytest.MonkeyPatch) -> None:
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
        def __enter__(self) -> "FakeClient":
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
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
            return minimal_args.output_csv

        return fetcher, writer

    def fake_run_pipeline(**kwargs: object) -> int:
        kwargs["stats_callback"]({"rows_total": 2, "rows_kept": 2, "rows_dropped": 0})
        return 0

    monkeypatch.setattr(get_assay_data.io, "read_ids", fake_read_ids)
    monkeypatch.setattr(get_assay_data, "ChemblClient", lambda *args, **kwargs: FakeClient())
    monkeypatch.setattr(get_assay_data, "ChunkFailureTracker", lambda: tracker)
    monkeypatch.setattr(get_assay_data.cl, "get_assays", lambda *args, **kwargs: pd.DataFrame({"assay_chembl_id": ["CHEMBL1"]}))
    monkeypatch.setattr(get_assay_data, "prepare_chunked_pipeline", fake_prepare_chunked_pipeline)
    monkeypatch.setattr(get_assay_data, "run_pipeline", fake_run_pipeline)

    exit_code = get_assay_data.run_chembl(cfg, minimal_args)

    assert exit_code == 0
    assert any(event == "assay_pipeline_done" for _, event, _ in logger_stub.events)
    if offset:
        assert any(event == "process_offset" for _, event, _ in logger_stub.events)
    assert any(event == "process_limit" for _, event, _ in logger_stub.events)


def test_run__skip_existing_returns_zero(
    cfg: Config,
    minimal_args: argparse.Namespace,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    minimal_args.skip_existing = True
    minimal_args.force = False
    minimal_args.output_csv.write_text("existing", encoding="utf-8")

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
        {"output": str(minimal_args.output_csv)},
    ) in logger_stub.events


def test_run__force_overrides_skip(
    cfg: Config,
    minimal_args: argparse.Namespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    minimal_args.skip_existing = True
    minimal_args.force = True
    minimal_args.output_csv.write_text("existing", encoding="utf-8")

    calls: list[str] = []

    def fake_run(cfg: Config, args: argparse.Namespace) -> int:
        calls.append("run")
        return 3

    monkeypatch.setattr(get_assay_data, "run_chembl", fake_run)

    exit_code = get_assay_data.run(cfg, minimal_args)

    assert exit_code == 3
    assert calls == ["run"]


def test_run__propagates_exit_code(cfg: Config, minimal_args: argparse.Namespace, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(get_assay_data, "run_chembl", lambda *_: 7)

    exit_code = get_assay_data.run(cfg, minimal_args)

    assert exit_code == 7


def test_build_parser__defaults() -> None:
    parser, _ = get_assay_data.build_parser()
    args = parser.parse_args([])

    assert args.input_csv == Path(get_assay_data.DEFAULT_INPUT_NAME)
    assert args.batch_size == parser.get_default("batch_size")
    assert callable(args.func)
