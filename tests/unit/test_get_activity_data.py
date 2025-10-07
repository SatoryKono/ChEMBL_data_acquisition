"""Unit tests for :mod:`scripts.get_activity_data`."""

from __future__ import annotations

import argparse
from collections.abc import Iterable
from pathlib import Path

import pandas as pd
import pytest

from library.pipelines.assay.chembl_assay import ACTIVITY_COLUMNS
from library.postprocessing import activity_extended
from scripts import get_activity_data


def _make_args(tmp_path: Path) -> argparse.Namespace:
    """Return an :class:`argparse.Namespace` populated with sane defaults."""

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("activity_id\nA1\n", encoding="utf-8")
    return argparse.Namespace(
        input_csv=input_csv,
        final_out=tmp_path / "output.csv",
        offset=0,
        workers=None,
        skip_existing=False,
        force=False,
        dry_run=False,
        invocation=None,
    )


class _DummyClient:
    """Minimal context manager emulating :class:`ChemblClient`."""

    def __init__(self, *args, **kwargs) -> None:  # pragma: no cover - interface compatibility
        pass

    def __enter__(self) -> _DummyClient:  # pragma: no cover - trivial helper
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:  # pragma: no cover - trivial helper
        return False

    def close(self) -> None:  # pragma: no cover - interface compatibility
        return None


class _RecordingLogger:
    """Capture structured log events emitted by :mod:`get_activity_data`."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def info(self, event: str, **kwargs: object) -> None:
        self.events.append(("info", event, dict(kwargs)))

    def warning(self, event: str, **kwargs: object) -> None:
        self.events.append(("warning", event, dict(kwargs)))

    def error(self, event: str, **kwargs: object) -> None:
        self.events.append(("error", event, dict(kwargs)))


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


def test_run_chembl__dry_run_short_circuits(cfg, tmp_path, monkeypatch) -> None:
    args = _make_args(tmp_path)
    cfg.activity.dry_run = True
    cfg.activity.limit = 7

    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 0
    assert ("info", "dry_run", {"limit": 7}) in [
        (level, event, context) for level, event, context in logger_stub.events
    ]
    completion_events = [
        event
        for _, event, _ in logger_stub.events
        if event.startswith("Completed get_activity_data pipeline:")
    ]
    assert completion_events
    assert "mode=dry_run" in completion_events[-1]


@pytest.mark.parametrize(
    "skip_existing, force, explicit_output, has_existing, expected_calls",
    [
        (True, False, False, True, 0),
        (True, True, False, True, 1),
        (False, False, True, True, 1),
        (False, False, True, False, 1),
    ],
)
def test_run__skip_existing_matrix(
    skip_existing: bool,
    force: bool,
    explicit_output: bool,
    has_existing: bool,
    expected_calls: int,
    cfg,
    tmp_path,
    monkeypatch,
) -> None:
    args = _make_args(tmp_path)
    args.skip_existing = skip_existing
    args.force = force

    if explicit_output:
        output_path = tmp_path / "explicit.csv"
        args.final_out = output_path
    else:
        args.final_out = None
        output_path = tmp_path / "default.csv"

        def fake_default_output_path(input_path: Path, _io_cfg) -> Path:
            assert input_path == args.input_csv
            return output_path

        monkeypatch.setattr(get_activity_data.io, "default_output_path", fake_default_output_path)

    if has_existing:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("activity_id\nA1\n", encoding="utf-8")

    call_counter: list[str] = []
    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr(
        get_activity_data,
        "run_chembl",
        lambda *_args, **_kwargs: call_counter.append("run") or 0,
    )

    exit_code = get_activity_data.run(cfg, args)

    assert exit_code == 0
    assert len(call_counter) == expected_calls
    summary_events = [
        event
        for _, event, _ in logger_stub.events
        if event.startswith("Completed get_activity_data pipeline:")
    ]
    if skip_existing and has_existing and not force:
        assert summary_events
        assert "mode=skip_existing" in summary_events[-1]
    else:
        assert not summary_events


@pytest.mark.parametrize(
    ("argv", "expected_substring"),
    [
        (["--limit", "-1"], "--limit must be zero or a positive integer"),
        (["--offset", "-2"], "--offset must be zero or a positive integer"),
        (["--workers", "0"], "chunk size must be a positive integer"),
    ],
)
def test_main__invalid_cli_options(argv, expected_substring, monkeypatch, tmp_path, capsys) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    full_args = ["--input", str(input_csv), *argv]

    monkeypatch.setattr(
        get_activity_data,
        "run",
        lambda *args, **kwargs: pytest.fail("run should not be invoked for invalid CLI"),
    )

    with pytest.raises(SystemExit) as exc_info:
        get_activity_data.main(full_args)

    assert exc_info.value.code != 0
    captured = capsys.readouterr()
    assert expected_substring in captured.err


@pytest.mark.parametrize(
    "argv", [["--limit"], ["--workers"], ["--timeout"]]
)
def test_main__missing_cli_option_values(argv, monkeypatch, tmp_path, capsys) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    args = ["--input", str(input_csv), *argv]

    monkeypatch.setattr(
        get_activity_data,
        "run",
        lambda *_, **__: pytest.fail("run should not execute when CLI value missing"),
    )

    with pytest.raises(SystemExit) as exc_info:
        get_activity_data.main(args)

    assert exc_info.value.code != 0
    stderr = capsys.readouterr().err
    option = argv[0]
    assert f"argument {option}: expected one argument" in stderr


def test_run_chembl__offset_and_workers(monkeypatch, cfg, tmp_path) -> None:
    args = _make_args(tmp_path)
    args.offset = 1
    cfg.activity.dry_run = False
    cfg.activity.limit = None
    cfg.activity.batch_size = 3
    cfg.activity.workers = 4
    cfg.system.doc_quality.enable = False

    monkeypatch.setattr(
        get_activity_data.io,
        "read_ids",
        lambda *_args, **_kwargs: iter(["ACT0", "ACT1", "ACT2"]),
    )

    captured: dict[str, int] = {}

    def fake_prepare_chunked_pipeline(*, fetch_config, fetch_chunk, csv_writer):
        captured["workers"] = fetch_config.workers
        captured["chunk_size"] = fetch_config.chunk_size

        def _fetcher() -> Iterable[pd.DataFrame]:
            yield pd.DataFrame(
                [
                    {
                        "activity_id": "ACT2",
                        "molecule_chembl_id": "CHEMBL2",
                        "assay_chembl_id": "ASSAY2",
                        "standard_value": 4.0,
                    }
                ]
            )

        def _writer(chunks: Iterable[pd.DataFrame], destination: Path, col_order, key_cols):
            frames = list(chunks)
            result = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
            destination.parent.mkdir(parents=True, exist_ok=True)
            result.to_csv(destination, index=False)
            return destination

        return _fetcher, _writer

    monkeypatch.setattr(get_activity_data, "prepare_chunked_pipeline", fake_prepare_chunked_pipeline)
    monkeypatch.setattr("library.orchestration.context.ChemblClient", _DummyClient)
    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)

    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 0
    assert captured["workers"] == max(1, cfg.activity.workers)
    events = [event for _, event, _ in logger_stub.events]
    assert "process_offset" in events


@pytest.mark.parametrize(
    "error_factory",
    [
        lambda: FileNotFoundError("missing input"),
        lambda: ValueError("malformed CSV"),
    ],
)
def test_run_chembl__read_ids_failures(error_factory, cfg, tmp_path, monkeypatch) -> None:
    args = _make_args(tmp_path)

    def _raise(*_args, **_kwargs):
        raise error_factory()

    monkeypatch.setattr(get_activity_data.io, "read_ids", _raise)
    monkeypatch.setattr("library.orchestration.context.ChemblClient", _DummyClient)
    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)

    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 1
    events = [event for _, event, _ in logger_stub.events]
    assert "read_fail" in events


def test_run_chembl__pipeline_failure_logs_error(cfg, tmp_path, monkeypatch) -> None:
    args = _make_args(tmp_path)

    monkeypatch.setattr(
        get_activity_data.io,
        "read_ids",
        lambda *_args, **_kwargs: iter(["ACT1"]),
    )
    monkeypatch.setattr("library.orchestration.context.ChemblClient", _DummyClient)

    def fake_prepare_chunked_pipeline(*, fetch_config, fetch_chunk, csv_writer):
        def _fetcher() -> Iterable[pd.DataFrame]:
            yield pd.DataFrame(
                [
                    {
                        "activity_id": "ACT1",
                        "molecule_chembl_id": "CHEMBL1",
                        "assay_chembl_id": "ASSAY1",
                        "standard_value": 1.0,
                    }
                ]
            )

        def _writer(
            chunks: Iterable[pd.DataFrame],
            destination: Path,
            col_order,
            key_cols,
        ) -> Path:
            dest = Path(destination)
            dest.parent.mkdir(parents=True, exist_ok=True)
            return dest

        return _fetcher, _writer

    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr(
        get_activity_data, "prepare_chunked_pipeline", fake_prepare_chunked_pipeline
    )

    def fake_run_pipeline(*, definition, fetcher, output_path, failure_path, **kwargs):
        del definition, output_path, failure_path, kwargs
        list(fetcher())
        return 1

    monkeypatch.setattr(get_activity_data, "run_pipeline", fake_run_pipeline)

    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 1
    error_events = [event for level, event, _ in logger_stub.events if level == "error"]
    assert "activity_pipeline_failed" in error_events


def test_main__dry_run_skip_limit(monkeypatch, tmp_path, capsys) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")

    monkeypatch.setattr(
        get_activity_data,
        "run",
        lambda *args, **kwargs: pytest.fail("run should not execute when limit=0"),
    )
    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)

    exit_code = get_activity_data.main(["--input", str(input_csv), "--limit", "0"])

    assert exit_code == 0
    events = [event for _, event, _ in logger_stub.events]
    assert "pipeline_skip_limit" in events


def test_activity_columns__cover_extended_requirements() -> None:
    required = activity_extended._REQUIRED_COLUMNS
    available = set(ACTIVITY_COLUMNS)
    assert required <= available


def test_ensure_extended_activity_columns__adds_defaults() -> None:
    frame = pd.DataFrame(
        {
            "activity_id": ["A1"],
            "molecule_chembl_id": ["M1"],
            "pchembl_value": [5.3],
            "bao_label": ["BAO:000001"],
        }
    )

    enriched = get_activity_data._ensure_extended_activity_columns(frame)

    assert enriched["activity_chembl_id"].tolist() == ["A1"]
    assert enriched.loc[0, "log_value"] == pytest.approx(5.3)
    assert enriched["compound_name"].isna().all()
    assert enriched["salt_chembl_id"].isna().all()


def test_ensure_extended_activity_columns__preserves_salt_ids() -> None:
    frame = pd.DataFrame(
        {
            "activity_id": ["A1", "A2"],
            "salt_chembl_id": ["CHEMBL_SALT1", pd.NA],
            "molecule_chembl_id": ["CHEMBL123", "CHEMBL456"],
        }
    )

    enriched = get_activity_data._ensure_extended_activity_columns(frame)

    assert enriched["salt_chembl_id"].tolist() == ["CHEMBL_SALT1", pd.NA]


def test_ensure_src_assay_id__fills_from_lookup() -> None:
    frame = pd.DataFrame(
        {
            "assay_chembl_id": ["ASSAY1", "ASSAY2", "ASSAY3"],
            "src_assay_id": ["", pd.NA, "custom"],
        }
    )

    enriched = get_activity_data._ensure_src_assay_id(
        frame,
        lookup={
            "ASSAY1": "SRC-1",
            "ASSAY2": "SRC-2",
        },
    )

    assert enriched["src_assay_id"].tolist() == ["SRC-1", "SRC-2", "custom"]
    assert pd.api.types.is_string_dtype(enriched["src_assay_id"])  # type: ignore[arg-type]


def test_ensure_src_assay_id__adds_column_when_missing() -> None:
    frame = pd.DataFrame({"assay_chembl_id": ["ASSAY1"]})

    enriched = get_activity_data._ensure_src_assay_id(frame, lookup={})

    assert "src_assay_id" in enriched.columns
    assert enriched["src_assay_id"].isna().all()
    assert pd.api.types.is_string_dtype(enriched["src_assay_id"])  # type: ignore[arg-type]
