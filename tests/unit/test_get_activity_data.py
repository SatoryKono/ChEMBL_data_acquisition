"""Unit tests for :mod:`scripts.get_activity_data`."""

from __future__ import annotations

import argparse
import importlib
import sys
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Condition, Event, Lock
from types import ModuleType, SimpleNamespace
from typing import Any

import pandas as pd
import pytest

from library.pipelines.assay.chembl_assay import ACTIVITY_COLUMNS
from library.pipelines.common import PipelineRunResult
from library.postprocessing import activity_extended
from scripts import get_activity_data


def test_command_module_run_chembl__delegates_to_entrypoint(monkeypatch):
    module_name = "library.cli.entrypoints.activity"
    stub_module = ModuleType(module_name)
    calls: dict[str, object] = {}

    def _run(cfg, args):
        calls["cfg"] = cfg
        calls["args"] = args
        return 123

    def _emit(*args, **kwargs):
        calls["emit"] = (args, kwargs)

    stub_module.run_chembl = _run  # type: ignore[attr-defined]
    stub_module._emit_completion_message = _emit  # type: ignore[attr-defined]
    stub_module.__all__ = ("run_chembl", "_emit_completion_message")

    import library.cli.activity_api as activity_api

    monkeypatch.setattr(activity_api, "_ENTRYPOINT_MODULE_CACHE", None)
    monkeypatch.setitem(sys.modules, module_name, stub_module)

    command_module = importlib.reload(
        importlib.import_module("library.cli.commands.get_activity_data")
    )

    result = command_module.run_chembl("cfg", "args")
    command_module._emit_completion_message(1, mode="test")

    assert result == 123
    assert calls["cfg"] == "cfg"
    assert calls["args"] == "args"
    assert calls["emit"] == ((1,), {"mode": "test"})


def test_wrapper_module__reflects_command_updates(monkeypatch):
    command_module = importlib.import_module("library.cli.commands.get_activity_data")
    wrapper_module = importlib.import_module("scripts.get_activity_data")

    assert wrapper_module is command_module

    sentinel = object()
    monkeypatch.setattr(command_module, "_sentinel_for_test", sentinel, raising=False)

    assert wrapper_module._sentinel_for_test is sentinel


def _make_pref_name_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "molecule_chembl_id": pd.Series(["CHEMBL1", "CHEMBL2"], dtype="string"),
            "molecule_pref_name": pd.Series([pd.NA, pd.NA], dtype="string"),
        }
    )


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


def test_build_parser__captures_run_id_env(monkeypatch, tmp_path):
    monkeypatch.setenv("CHEMBL_DA_RUN_ID", "env-run")
    parser, log_cfg = get_activity_data.build_parser()

    args = parser.parse_args(
        ["--input", str(tmp_path / "activity.csv"), "--limit", "0"]
    )

    assert args.run_id == "env-run"
    assert log_cfg.run_id == "env-run"


class _DummyClient:
    """Minimal context manager emulating :class:`ChemblClient`."""

    def __init__(
        self, *args, **kwargs
    ) -> None:  # pragma: no cover - interface compatibility
        pass

    def __enter__(self) -> _DummyClient:  # pragma: no cover - trivial helper
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:  # pragma: no cover - trivial helper
        return False

    def close(self) -> None:  # pragma: no cover - interface compatibility
        return None


class _ApiStub:
    """Lightweight stand-in for :class:`ApiCfg` supporting ``model_copy``."""

    def __init__(self, **attributes: Any) -> None:
        self.__dict__.update(attributes)

    def model_copy(
        self,
        *,
        update: Mapping[str, Any] | None = None,
        deep: bool | None = None,  # noqa: ARG002 - signature compatibility
    ) -> _ApiStub:
        data = dict(self.__dict__)
        if update:
            data.update(update)
        return _ApiStub(**data)


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
def test_args_invocation__cases(
    invocation: Iterable[object] | None, expected: tuple[str, ...]
) -> None:
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
        payload
        for _, event, payload in logger_stub.events
        if event == "activity_pipeline_completion"
    ]
    assert completion_events
    assert completion_events[-1]["mode"] == "dry_run"


def test_ensure_molecule_pref_name__concurrent_single_fetch(monkeypatch) -> None:
    cfg = SimpleNamespace(
        api=_ApiStub(retries=1, backoff_factor=0.1),
        testitem=SimpleNamespace(
            fields=["pref_name"],
            batch_size=10,
            timeout=5.0,
            request_limit=5,
            retries=None,
            backoff_factor=None,
        ),
    )

    frame = _make_pref_name_frame()
    cache: dict[str, str | None] = {}
    cache_lock = Lock()
    cache_condition = Condition(cache_lock)

    call_records: list[tuple[str, ...]] = []
    ready = Event()
    proceed = Event()

    def fake_get_testitem(
        identifiers: Sequence[str],
        *,
        cfg,
        client,
        chunk_size,
        timeout,
        fields,
        page_limit,
    ) -> pd.DataFrame:
        call_records.append(tuple(identifiers))
        ready.set()
        if not proceed.wait(timeout=1):
            pytest.fail("proceed event was not set")
        return pd.DataFrame(
            {
                "molecule_chembl_id": pd.Series(
                    [str(identifier) for identifier in identifiers], dtype="string"
                ),
                "pref_name": pd.Series(
                    [f"name-{identifier}" for identifier in identifiers], dtype="string"
                ),
            }
        )

    monkeypatch.setattr(get_activity_data.cl, "get_testitem", fake_get_testitem)

    def worker() -> pd.DataFrame:
        return get_activity_data._ensure_molecule_pref_name(
            frame,
            cfg=cfg,
            client=object(),
            cache=cache,
            cache_condition=cache_condition,
        )

    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(worker) for _ in range(4)]
        if not ready.wait(timeout=1):
            proceed.set()
            pytest.fail("cl.get_testitem was not invoked")
        proceed.set()
        results = [future.result() for future in futures]

    expected_series = pd.Series(["name-CHEMBL1", "name-CHEMBL2"], dtype="string")
    for result in results:
        pd.testing.assert_series_equal(
            result["molecule_pref_name"],
            expected_series.reindex(result.index),
            check_names=False,
        )

    identifier_calls = Counter()
    for entry in call_records:
        identifier_calls.update(entry)

    assert identifier_calls == Counter({"CHEMBL1": 1, "CHEMBL2": 1})
    assert len(call_records) == 1

    with cache_lock:
        assert cache["CHEMBL1"] == "name-CHEMBL1"
        assert cache["CHEMBL2"] == "name-CHEMBL2"


def test_ensure_molecule_pref_name__applies_testitem_retry_overrides(
    monkeypatch,
) -> None:
    cfg = SimpleNamespace(
        api=_ApiStub(retries=2, backoff_factor=0.5),
        testitem=SimpleNamespace(
            fields=["pref_name"],
            batch_size=2,
            timeout=5.0,
            request_limit=5,
            retries=7,
            backoff_factor=3.5,
        ),
    )

    frame = _make_pref_name_frame()
    cache: dict[str, str | None] = {}
    cache_lock = Lock()
    cache_condition = Condition(cache_lock)
    captured_cfg: dict[str, Any] = {}

    def fake_get_testitem(
        identifiers: Sequence[str],
        *,
        cfg,
        client,
        chunk_size,
        timeout,
        fields,
        page_limit,
    ) -> pd.DataFrame:
        captured_cfg["cfg"] = cfg
        return pd.DataFrame(
            {
                "molecule_chembl_id": pd.Series(
                    [str(identifier) for identifier in identifiers], dtype="string"
                ),
                "pref_name": pd.Series(
                    [f"name-{identifier}" for identifier in identifiers], dtype="string"
                ),
            }
        )

    monkeypatch.setattr(get_activity_data.cl, "get_testitem", fake_get_testitem)

    enriched = get_activity_data._ensure_molecule_pref_name(
        frame,
        cfg=cfg,
        client=object(),
        cache=cache,
        cache_condition=cache_condition,
    )

    assert "cfg" in captured_cfg
    effective_cfg = captured_cfg["cfg"]
    assert isinstance(effective_cfg, _ApiStub)
    assert effective_cfg is not cfg.api
    assert effective_cfg.retries == cfg.testitem.retries
    assert effective_cfg.backoff_factor == cfg.testitem.backoff_factor

    expected_series = pd.Series(["name-CHEMBL1", "name-CHEMBL2"], dtype="string")
    pd.testing.assert_series_equal(
        enriched["molecule_pref_name"], expected_series, check_names=False
    )


def test_ensure_molecule_pref_name__requests_minimal_field_set(monkeypatch) -> None:
    cfg = SimpleNamespace(
        api=_ApiStub(),
        testitem=SimpleNamespace(
            fields=["molecule_chembl_id", "pref_name", "molecule_properties"],
            batch_size=2,
            timeout=5.0,
            request_limit=5,
            retries=None,
            backoff_factor=None,
        ),
    )

    frame = _make_pref_name_frame()
    cache: dict[str, str | None] = {}
    cache_lock = Lock()
    cache_condition = Condition(cache_lock)
    captured_fields: list[Sequence[str]] = []

    def fake_get_testitem(
        identifiers: Sequence[str],
        *,
        cfg,
        client,
        chunk_size,
        timeout,
        fields,
        page_limit,
    ) -> pd.DataFrame:
        captured_fields.append(tuple(fields))
        return pd.DataFrame(
            {
                "molecule_chembl_id": pd.Series(
                    [str(identifier) for identifier in identifiers], dtype="string"
                ),
                "pref_name": pd.Series(
                    [f"name-{identifier}" for identifier in identifiers], dtype="string"
                ),
            }
        )

    monkeypatch.setattr(get_activity_data.cl, "get_testitem", fake_get_testitem)

    enriched = get_activity_data._ensure_molecule_pref_name(
        frame,
        cfg=cfg,
        client=object(),
        cache=cache,
        cache_condition=cache_condition,
    )

    assert captured_fields == [
        ("molecule_chembl_id", "molecule_properties", "pref_name")
    ]

    expected_series = pd.Series(["name-CHEMBL1", "name-CHEMBL2"], dtype="string")
    pd.testing.assert_series_equal(
        enriched["molecule_pref_name"], expected_series, check_names=False
    )


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

        monkeypatch.setattr(
            get_activity_data.io, "default_output_path", fake_default_output_path
        )

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
        payload
        for _, event, payload in logger_stub.events
        if event == "activity_pipeline_completion"
    ]
    if skip_existing and has_existing and not force:
        assert summary_events
        assert summary_events[-1]["mode"] == "skip_existing"
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
def test_main__invalid_cli_options(
    argv, expected_substring, monkeypatch, tmp_path, capsys
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    full_args = ["--input", str(input_csv), *argv]

    monkeypatch.setattr(
        get_activity_data,
        "run",
        lambda *args, **kwargs: pytest.fail(
            "run should not be invoked for invalid CLI"
        ),
    )

    with pytest.raises(SystemExit) as exc_info:
        get_activity_data.main(full_args)

    assert exc_info.value.code != 0
    captured = capsys.readouterr()
    assert expected_substring in captured.err


@pytest.mark.parametrize("argv", [["--limit"], ["--workers"], ["--timeout"]])
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

    def fake_run_pipeline(**kwargs):
        fetch_config = kwargs["fetch_config"]
        captured["workers"] = fetch_config.workers
        captured["chunk_size"] = fetch_config.chunk_size
        output_path = kwargs["output_path"]
        output_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            [
                {
                    "activity_id": "ACT2",
                    "molecule_chembl_id": "CHEMBL2",
                    "assay_chembl_id": "ASSAY2",
                    "standard_value": 4.0,
                }
            ]
        ).to_csv(output_path, index=False)
        return PipelineRunResult(exit_code=0, output_path=output_path, written=True)

    monkeypatch.setattr(
        "library.pipelines.activity.run.run_activity_pipeline",
        fake_run_pipeline,
    )
    monkeypatch.setattr("library.orchestration.context.ChemblClient", _DummyClient)
    monkeypatch.setattr(
        get_activity_data,
        "process_activity_extended",
        lambda *, input_path, **__: input_path,
    )
    monkeypatch.setattr(
        get_activity_data,
        "_generate_activity_postprocess_metrics",
        lambda *_, **__: (None, None),
    )
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
def test_run_chembl__read_ids_failures(
    error_factory, cfg, tmp_path, monkeypatch
) -> None:
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

    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr(
        "library.pipelines.activity.run.run_activity_pipeline",
        lambda **kwargs: PipelineRunResult(
            exit_code=1, output_path=kwargs["output_path"], written=None
        ),
    )

    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 1
    error_events = [event for level, event, _ in logger_stub.events if level == "error"]
    assert "activity_pipeline_failed" in error_events


def test_prepare_activity_context__skip_read_avoids_io(
    cfg, tmp_path, monkeypatch
) -> None:
    args = _make_args(tmp_path)
    cfg.activity.limit = None

    monkeypatch.setattr(
        get_activity_data.io,
        "read_ids",
        lambda *_args, **_kwargs: pytest.fail(
            "read_ids should not execute when skip_read=True"
        ),
    )

    context = get_activity_data.prepare_activity_context(cfg, args, skip_read=True)

    assert context is not None
    assert context.limit is None
    assert list(context.limited_ids) == []
    assert context.processed_ids == 0


def test_prepare_activity_context__limit_and_offset(cfg, tmp_path, monkeypatch) -> None:
    args = _make_args(tmp_path)
    args.offset = 1
    cfg.activity.limit = 2

    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr(
        get_activity_data.io,
        "read_ids",
        lambda *_args, **_kwargs: iter(["ACT0", "ACT1", "ACT2", "ACT3"]),
    )

    context = get_activity_data.prepare_activity_context(cfg, args)

    assert context is not None
    assert list(context.limited_ids) == ["ACT1", "ACT2"]
    assert context.processed_ids == 2
    events = [event for _, event, _ in logger_stub.events]
    assert "process_offset" in events


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


def test_main__missing_input_fails_fast(tmp_path, capsys) -> None:
    missing = tmp_path / "absent.csv"

    exit_code = get_activity_data.main(["--input", str(missing)])

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "does not exist" in captured.err
    assert "Provide --input" in captured.err


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


@pytest.mark.unit
def test_run_chembl__delegates_to_run_activity_pipeline(monkeypatch):
    module = importlib.import_module("library.cli.commands.get_activity_data")

    called = {}

    def fake_pipeline(*args, **kwargs):
        called["args"] = args
        called["kwargs"] = kwargs
        return "sentinel"

    monkeypatch.setattr(module, "run_activity_pipeline", fake_pipeline)

    result = module.run_chembl("cfg", option=True)

    assert result == "sentinel"
    assert called["args"] == ("cfg",)
    assert called["kwargs"] == {"option": True}


@pytest.mark.unit
def test_run_chembl__exposed_via_all():
    module = importlib.import_module("library.cli.commands.get_activity_data")

    assert "run_chembl" in module.__all__
    assert module.run_chembl is not None
