"""Unit tests for :mod:`scripts.get_activity_data`."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd
import pytest

from scripts import get_activity_data


def _make_args(tmp_path: Path) -> argparse.Namespace:
    """Return an :class:`argparse.Namespace` populated with sane defaults."""

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("activity_id\nA1\n", encoding="utf-8")
    return argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "output.csv",
        offset=0,
        workers=None,
        skip_existing=False,
        force=False,
        dry_run=False,
        invocation=None,
    )


@pytest.mark.parametrize(
    "invocation, expected",
    [
        (None, (get_activity_data.PROGRAM_NAME,)),
        ((), ()),
        ([], ()),
        (["run"], ("run",)),
        (("fetch", 2), ("fetch", "2")),
        ([Path("/tmp/file.csv")], ("/tmp/file.csv",)),
        ([""], ("",)),
        (["--limit", 5], ("--limit", "5")),
        (["Σ", Path("relative")], ("Σ", "relative")),
        ([42.0, 0], ("42.0", "0")),
    ],
)
def test_args_invocation__cases(invocation: Iterable[object] | None, expected: tuple[str, ...]) -> None:
    namespace = argparse.Namespace(invocation=invocation)
    assert get_activity_data._args_invocation(namespace) == expected


def test_run_chembl__invalid_limit_logs_error(cfg, tmp_path, caplog) -> None:
    args = _make_args(tmp_path)
    cfg.activity.limit = -5

    caplog.set_level("ERROR")
    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 1


def test_run_chembl__dry_run_short_circuits(cfg, tmp_path, caplog) -> None:
    args = _make_args(tmp_path)
    cfg.activity.dry_run = True
    cfg.activity.limit = 7

    caplog.set_level("INFO")
    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 0


def test_run__skip_existing_respects_flags(cfg, tmp_path, monkeypatch) -> None:
    args = _make_args(tmp_path)
    args.output_csv = None
    args.skip_existing = True
    args.force = False

    calls: list[str] = []

    def fake_default_output_path(input_path: Path, _io_cfg) -> Path:  # pragma: no cover - exercised via tests
        assert input_path == args.input_csv
        destination = tmp_path / "output.csv"
        destination.write_text("existing", encoding="utf-8")
        return destination

    monkeypatch.setattr(get_activity_data.io, "default_output_path", fake_default_output_path)
    monkeypatch.setattr(get_activity_data, "run_chembl", lambda *_: calls.append("run") or 0)

    exit_code = get_activity_data.run(cfg, args)

    assert exit_code == 0
    assert calls == []


class _MemoryLogger:
    """Capture structured log events emitted by the activity pipeline."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def info(self, event: str, **payload: object) -> None:
        self.events.append(("info", event, dict(payload)))

    def warning(self, event: str, **payload: object) -> None:
        self.events.append(("warning", event, dict(payload)))

    def error(self, event: str, **payload: object) -> None:
        self.events.append(("error", event, dict(payload)))

    def debug(self, event: str, **payload: object) -> None:  # pragma: no cover - rarely used
        self.events.append(("debug", event, dict(payload)))

    def bind(self, **_: object) -> "_MemoryLogger":  # pragma: no cover - logging compat
        return self


@pytest.fixture()
def logger_stub(monkeypatch: pytest.MonkeyPatch) -> _MemoryLogger:
    logger = _MemoryLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger)
    return logger


class _DummyClient:
    def __enter__(self) -> "_DummyClient":  # pragma: no cover - simple context helper
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # pragma: no cover - no special cleanup
        return None


class _DummyTracker:
    def __init__(self) -> None:
        self.stats: dict[str, object] = {"failures": 0}
        self.saved_path: Path | None = None

    def add_failure(self, *_: object, **__: object) -> None:
        return None

    def save(self, path: Path, *, cfg: object) -> None:
        self.saved_path = path


def _setup_pipeline(
    monkeypatch: pytest.MonkeyPatch,
    cfg,
    args: argparse.Namespace,
    *,
    exit_code: int = 0,
    stats: dict[str, int] | None = None,
) -> tuple[dict[str, object], _DummyTracker]:
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        get_activity_data.io,
        "read_ids",
        lambda path, *, column, cfg: iter(["CHEMBL1", "CHEMBL2", "CHEMBL3"]),
    )
    monkeypatch.setattr(get_activity_data, "ChemblClient", lambda *a, **k: _DummyClient())
    tracker = _DummyTracker()
    monkeypatch.setattr(get_activity_data, "ChunkFailureTracker", lambda: tracker)
    monkeypatch.setattr(get_activity_data, "configure_activity_schema", lambda *_: None)
    monkeypatch.setattr(
        get_activity_data,
        "compute_activity_bounds",
        lambda frame, *_: frame.assign(bounds_applied=True),
    )
    monkeypatch.setattr(
        get_activity_data,
        "apply_activity_annotations",
        lambda frame, **__: frame.assign(annotations_applied=True),
    )
    monkeypatch.setattr(
        get_activity_data,
        "write_csv_chunks_deterministic",
        lambda chunks, destination, **__: destination,
    )

    def fake_prepare_chunked_pipeline(*, fetch_config, fetch_chunk, csv_writer):
        ids = list(fetch_config.ids)
        captured["fetch_workers"] = fetch_config.workers
        captured["fetch_chunk_size"] = fetch_config.chunk_size
        captured["consumed_ids"] = tuple(ids)

        frames = [pd.DataFrame({"activity_id": ids})]

        def fetcher() -> Iterable[pd.DataFrame]:
            yield from frames

        def writer(
            chunks: Iterable[pd.DataFrame],
            destination: Path,
            col_order: Iterable[str],
            key_cols: Iterable[str],
        ) -> Path:
            captured["writer_destination"] = destination
            captured["writer_columns"] = tuple(col_order)
            captured["writer_key_cols"] = tuple(key_cols)
            captured["written_rows"] = sum(len(chunk) for chunk in chunks)
            return destination

        return fetcher, writer

    monkeypatch.setattr(get_activity_data, "prepare_chunked_pipeline", fake_prepare_chunked_pipeline)
    monkeypatch.setattr(
        get_activity_data.cl,
        "get_activities",
        lambda chunk_ids, **__: pd.DataFrame({"activity_id": list(chunk_ids)}),
    )

    def fake_run_pipeline(**kwargs: object) -> int:
        fetcher = kwargs["fetcher"]
        metadata_hooks = list(kwargs["metadata_hooks"])
        writer = kwargs["writer"]
        output_path = kwargs["output_path"]
        key_columns = kwargs["key_columns"]
        schema = kwargs["schema"]
        stats_callback = kwargs["stats_callback"]

        frames_after_hooks: list[pd.DataFrame] = []
        for frame in fetcher():
            for hook in metadata_hooks:
                frame = hook(frame)
            frames_after_hooks.append(frame)
        captured["frames_after_hooks"] = frames_after_hooks

        final_stats = stats or {
            "rows_total": sum(len(frame) for frame in frames_after_hooks),
            "rows_kept": sum(len(frame) for frame in frames_after_hooks),
            "rows_dropped": 0,
        }
        stats_callback(final_stats)

        writer(
            frames_after_hooks,
            output_path,
            list(getattr(schema, "columns").keys()),
            key_columns,
        )
        return exit_code

    monkeypatch.setattr(get_activity_data, "run_pipeline", fake_run_pipeline)

    return captured, tracker


def test_run_chembl__read_ids_value_error_logs(cfg, tmp_path, monkeypatch, logger_stub) -> None:
    args = _make_args(tmp_path)

    def fake_read_ids(path: Path, *, column: str, cfg) -> Iterable[str]:
        raise ValueError("bad csv")

    monkeypatch.setattr(get_activity_data.io, "read_ids", fake_read_ids)

    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 1
    assert (
        "error",
        "read_fail",
        {"error": "bad csv", "path": str(args.input_csv)},
    ) in logger_stub.events


def test_run_chembl__offset_limit_and_success_logging(
    cfg, tmp_path, monkeypatch, logger_stub
) -> None:
    args = _make_args(tmp_path)
    args.offset = 1
    args.workers = 2
    cfg.activity.limit = 2

    captured, tracker = _setup_pipeline(monkeypatch, cfg, args)

    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 0
    assert captured["consumed_ids"] == ("CHEMBL2", "CHEMBL3")
    assert tracker.saved_path == args.output_csv.with_name(
        f"{args.output_csv.stem}_fetch_failures.csv"
    )
    assert captured["written_rows"] == 2
    assert captured["fetch_workers"] == max(1, cfg.activity.workers)

    events = logger_stub.events
    assert ("info", "process_offset", {"offset": 1}) in events
    assert ("info", "process_limit", {"limit": 2}) in events
    assert any(event == "records_dropped" for _, event, _ in events)
    assert any(event == "activity_pipeline_done" for _, event, _ in events)


def test_run_chembl__pipeline_failure_emits_error(
    cfg, tmp_path, monkeypatch, logger_stub
) -> None:
    args = _make_args(tmp_path)
    cfg.activity.limit = None

    captured, tracker = _setup_pipeline(
        monkeypatch,
        cfg,
        args,
        exit_code=5,
        stats={"rows_total": 3, "rows_kept": 2, "rows_dropped": 1},
    )

    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 5
    assert tracker.saved_path == args.output_csv.with_name(
        f"{args.output_csv.stem}_fetch_failures.csv"
    )

    events = logger_stub.events
    failure_events = [payload for level, name, payload in events if name == "activity_pipeline_failed"]
    assert failure_events
    assert failure_events[0]["exit_code"] == 5
    assert failure_events[0]["processed"] == len(captured["consumed_ids"])


def test_run_chembl__worker_count_floors_to_one(cfg, tmp_path, monkeypatch, logger_stub) -> None:
    args = _make_args(tmp_path)
    cfg.activity.limit = 1
    cfg.activity.workers = 0

    captured, _tracker = _setup_pipeline(monkeypatch, cfg, args)

    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 0
    assert captured["fetch_workers"] == 1


def test_run_chembl__metadata_hooks_fill_required_columns(
    cfg, tmp_path, monkeypatch, logger_stub
) -> None:
    args = _make_args(tmp_path)
    cfg.activity.limit = 1

    captured, _ = _setup_pipeline(monkeypatch, cfg, args)

    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 0
    frames = captured.get("frames_after_hooks", [])
    assert frames
    frame = frames[0]
    for column in get_activity_data._ACTIVITY_REQUIRED_COLUMNS:
        assert column in frame.columns
