"""Integration tests for the activity pipeline CLI helpers."""

from __future__ import annotations

import argparse
import json
import warnings
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import requests
import yaml

from config.paths import DICTIONARY_DIR
from library.cli.commands import get_activity_data as command_activity
from library.config import Config
from library.resources.dictionaries import get_resource
from scripts import get_activity_data


class _DummyChemblClient:
    """Minimal stand-in for :class:`ChemblClient` used in tests."""

    def __init__(
        self, *args, **kwargs
    ) -> None:  # pragma: no cover - interface compatibility
        pass

    def __enter__(self) -> _DummyChemblClient:  # pragma: no cover - trivial
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:  # pragma: no cover - trivial
        return False

    def close(self) -> None:  # pragma: no cover - interface compatibility
        return None


class _RecordingLogger:
    """Capture structured events emitted during the pipeline run."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def _record(
        self,
        level: str,
        event: str,
        args: tuple[object, ...],
        kwargs: dict[str, object],
    ) -> None:
        payload = dict(kwargs)
        if args:
            payload["args"] = args
        self.events.append((level, event, payload))

    def debug(self, event: str, *args: object, **kwargs: object) -> None:
        self._record("debug", event, args, dict(kwargs))

    def info(self, event: str, *args: object, **kwargs: object) -> None:
        self._record("info", event, args, dict(kwargs))

    def warning(self, event: str, *args: object, **kwargs: object) -> None:
        self._record("warning", event, args, dict(kwargs))

    def error(self, event: str, *args: object, **kwargs: object) -> None:
        self._record("error", event, args, dict(kwargs))


@pytest.fixture()
def activity_resource_dir(snapshot_resource: Path) -> Path:
    return snapshot_resource / "activity_pipeline"


def _copy_resource(resource_dir: Path, name: str, destination: Path) -> Path:
    target = destination / name
    target.write_text(
        (resource_dir / name).read_text(encoding="utf-8"), encoding="utf-8"
    )
    return target


@dataclass
class _FetchCapture:
    activities: list[tuple[str, ...]]
    testitems: list[tuple[str, ...]]


def _install_fetch_stubs(
    monkeypatch: pytest.MonkeyPatch,
    frame: pd.DataFrame,
    *,
    testitem_frame: pd.DataFrame | None = None,
    testitem_error: Callable[[], Exception] | Exception | None = None,
) -> _FetchCapture:
    captured_activities: list[tuple[str, ...]] = []
    captured_testitems: list[tuple[str, ...]] = []

    def _fake_get_activities(chunk_ids: Iterable[str], **_: object) -> pd.DataFrame:
        identifiers = [str(item) for item in chunk_ids]
        captured_activities.append(tuple(identifiers))
        mask = frame["activity_id"].astype(str).isin(identifiers)
        result = frame.loc[mask].copy().reset_index(drop=True)
        if identifiers:
            result = result.drop_duplicates(subset=["activity_id"], keep="first")
        return result

    monkeypatch.setattr(get_activity_data.cl, "get_activities", _fake_get_activities)
    if testitem_frame is None:
        testitem_frame = pd.DataFrame(columns=["molecule_chembl_id", "pref_name"])

    def _fake_get_testitem(chunk_ids: Iterable[str], **_: object) -> pd.DataFrame:
        identifiers = [str(item) for item in chunk_ids]
        captured_testitems.append(tuple(identifiers))
        if testitem_error is not None:
            exc = testitem_error() if callable(testitem_error) else testitem_error
            raise exc
        if identifiers:
            mask = testitem_frame["molecule_chembl_id"].astype(str).isin(identifiers)
            return testitem_frame.loc[mask].copy().reset_index(drop=True)
        return testitem_frame.iloc[0:0].copy()

    monkeypatch.setattr(get_activity_data.cl, "get_testitem", _fake_get_testitem)
    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        _DummyChemblClient,
    )
    return _FetchCapture(captured_activities, captured_testitems)


def _install_writer_stub(
    monkeypatch: pytest.MonkeyPatch,
) -> list[tuple[Path, pd.DataFrame]]:
    written: list[tuple[Path, pd.DataFrame]] = []

    def _writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        *,
        key_cols: Iterable[str] | None,
        col_order: Iterable[str] | None,
        chunksize: int,
        sort_chunksize: int | None = None,
        sep: str = ",",
        encoding: str = "utf-8",
        cfg=None,
        **_: object,
    ) -> Path:
        collected = [chunk.copy() for chunk in chunks]
        if collected:
            result = pd.concat(collected, ignore_index=True)
        else:
            result = pd.DataFrame(columns=list(col_order or []))
        if col_order:
            order = [str(col) for col in col_order if str(col) in result.columns]
            extras = sorted(col for col in result.columns if col not in order)
            result = result.reindex(columns=order + extras)
        destination_path = Path(destination)
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(destination_path, index=False, sep=sep, encoding=encoding)
        written.append((destination_path, result))
        return destination_path

    monkeypatch.setattr(get_activity_data, "write_csv_chunks_deterministic", _writer)
    return written


def _configure_cfg(cfg) -> None:
    cfg.activity.limit = None
    cfg.activity.dry_run = False
    cfg.activity.batch_size = 5
    cfg.activity.workers = 1
    cfg.system.doc_quality.enable = False
    cfg.activity_enrichment.action_type.enabled = False
    cfg.activity_enrichment.action_type.log_missing = False
    cfg.activity_enrichment.activity_properties.enabled = False


def _make_args(input_csv: Path, output_csv: Path) -> argparse.Namespace:
    return argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        offset=0,
        workers=None,
        skip_existing=False,
        force=False,
        dry_run=False,
        invocation=None,
    )


def _activity_options_from_args(
    args: argparse.Namespace,
) -> command_activity.ActivityCommandOptions:
    """Construct :class:`ActivityCommandOptions` mirroring CLI semantics."""

    return command_activity.ActivityCommandOptions(
        input_csv=args.input_csv,
        output_csv=args.output_csv,
        final_output=args.output_csv,
        limit=getattr(args, "limit", None),
        offset=getattr(args, "offset", 0),
        timeout=getattr(args, "timeout", None),
        batch_size=getattr(args, "batch_size", None),
        workers=getattr(args, "workers", None),
        dry_run=getattr(args, "dry_run", False),
        skip_existing=getattr(args, "skip_existing", False),
        force=getattr(args, "force", False),
        invocation=getattr(args, "invocation", None),
    )


def _patch_activity_loggers(
    monkeypatch: pytest.MonkeyPatch, logger_stub: _RecordingLogger
) -> None:
    """Ensure both CLI entry points emit logs to the stub."""

    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr(command_activity, "logger", logger_stub)
    monkeypatch.setattr("library.common.log.logger", logger_stub)
    monkeypatch.setattr("library.pipelines.activity.runner.logger", logger_stub)


def _invoke_activity_runner(
    cfg,
    args: argparse.Namespace,
    variant: str,
    monkeypatch: pytest.MonkeyPatch,
    *,
    runner: Callable[[Config, argparse.Namespace], int],
) -> int:
    """Execute the orchestrator via the requested variant."""

    if variant == "cli":
        monkeypatch.setattr(command_activity, "run_chembl", runner)
        return command_activity.run(cfg, args)

    helper_options = _activity_options_from_args(args)
    return command_activity.run_activity_pipeline(
        cfg,
        helper_options,
        runner=runner,
        emit_completion_message=get_activity_data._emit_completion_message,
    )


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
@pytest.mark.parametrize("runner_variant", ["cli", "api"])
def test_activity_pipeline__timeout_clamped_when_below_minimum(
    cfg, tmp_path, monkeypatch, runner_variant
) -> None:
    _configure_cfg(cfg)
    cfg.activity.timeout = get_activity_data.MIN_ACTIVITY_TIMEOUT - 5

    input_csv = tmp_path / "ids.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    output_csv = tmp_path / "activities.csv"

    logger_stub = _RecordingLogger()
    _patch_activity_loggers(monkeypatch, logger_stub)

    captured_timeout: dict[str, float] = {}

    def _run_chembl_stub(passed_cfg, _passed_args):
        captured_timeout["timeout"] = float(passed_cfg.activity.timeout)
        return 0

    args = _make_args(input_csv, output_csv)

    exit_code = _invoke_activity_runner(
        cfg,
        args,
        runner_variant,
        monkeypatch,
        runner=_run_chembl_stub,
    )

    warning_events = [
        event for level, event, _ in logger_stub.events if level == "warning"
    ]
    assert "activity_timeout_clamped" in warning_events
    assert captured_timeout["timeout"] == pytest.approx(
        get_activity_data.MIN_ACTIVITY_TIMEOUT
    )
    assert cfg.activity.timeout == pytest.approx(
        get_activity_data.MIN_ACTIVITY_TIMEOUT - 5
    )
    assert exit_code == 0


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
@pytest.mark.parametrize("runner_variant", ["cli", "api"])
def test_activity_pipeline__warns_when_retry_disabled(
    cfg, tmp_path, monkeypatch, runner_variant
) -> None:
    _configure_cfg(cfg)
    cfg.retry.max_attempts = 1

    input_csv = tmp_path / "ids.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    output_csv = tmp_path / "activities.csv"

    logger_stub = _RecordingLogger()
    _patch_activity_loggers(monkeypatch, logger_stub)

    def _run_stub(passed_cfg, _args):
        assert passed_cfg.retry.max_attempts == 1
        return 0

    args = _make_args(input_csv, output_csv)

    exit_code = _invoke_activity_runner(
        cfg,
        args,
        runner_variant,
        monkeypatch,
        runner=_run_stub,
    )

    warning_events = [
        event for level, event, _ in logger_stub.events if level == "warning"
    ]
    assert "activity_retry_disabled" in warning_events
    assert exit_code == 0


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
@pytest.mark.parametrize("runner_variant", ["cli", "api"])
def test_activity_pipeline__warns_when_api_retries_disabled(
    cfg, tmp_path, monkeypatch, runner_variant
) -> None:
    _configure_cfg(cfg)
    cfg.api.retries = 0

    input_csv = tmp_path / "ids.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    output_csv = tmp_path / "activities.csv"

    logger_stub = _RecordingLogger()
    _patch_activity_loggers(monkeypatch, logger_stub)

    def _run_stub(passed_cfg, _args):
        assert passed_cfg.api.retries == 0
        return 0

    args = _make_args(input_csv, output_csv)

    exit_code = _invoke_activity_runner(
        cfg,
        args,
        runner_variant,
        monkeypatch,
        runner=_run_stub,
    )

    warning_events = [
        event for level, event, _ in logger_stub.events if level == "warning"
    ]
    assert "activity_api_retry_disabled" in warning_events
    assert exit_code == 0


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_activity_pipeline__fallback_postprocess_report_created(
    activity_resource_dir: Path, cfg, tmp_path, monkeypatch
) -> None:
    _configure_cfg(cfg)
    cfg.io.output_dir = tmp_path
    cfg.io.csv_sep = ","
    cfg.io.csv_encoding = "utf-8"

    input_csv = _copy_resource(activity_resource_dir, "ids_happy.csv", tmp_path)
    output_csv = tmp_path / "activities.csv"
    chunk_df = pd.read_csv(activity_resource_dir / "chunk_happy.csv")

    monkeypatch.setattr(
        get_activity_data,
        "process_activity_extended",
        lambda **_: None,
    )

    captured = _install_fetch_stubs(monkeypatch, chunk_df)
    written = _install_writer_stub(monkeypatch)

    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr("library.validation.logger", logger_stub)

    monkeypatch.setattr(get_activity_data, "_emit_completion_message", lambda **_: None)

    class _StubMetrics:
        def __init__(self) -> None:
            self.pipeline_version = "fallback-version"
            self.validation = SimpleNamespace(schema="activity-schema")

        def summary(self) -> dict[str, object]:
            return {
                "rows": 3,
                "columns": len(chunk_df.columns),
                "duration_s": 0.0,
                "steps": 2,
            }

    stub_metrics = _StubMetrics()
    collect_calls: list[dict[str, object]] = []
    original_collect = get_activity_data.collect_postprocess_metrics

    def _fake_collect_postprocess_metrics(**kwargs):
        collect_calls.append(kwargs)
        return stub_metrics, None

    monkeypatch.setattr(
        get_activity_data,
        "collect_postprocess_metrics",
        _fake_collect_postprocess_metrics,
    )

    try:
        args = _make_args(input_csv, output_csv)
        exit_code = get_activity_data.run(cfg, args)
    finally:
        monkeypatch.setattr(
            get_activity_data,
            "collect_postprocess_metrics",
            original_collect,
        )

    assert exit_code == 0
    assert collect_calls
    call_kwargs = collect_calls[0]
    assert call_kwargs["table"] == "activity"
    assert Path(call_kwargs["output_path"]) == output_csv

    fallback_path = tmp_path / "activity.postprocess.report.json"
    assert fallback_path.exists()

    payload = json.loads(fallback_path.read_text(encoding="utf-8"))
    assert payload["table"] == "activity"
    assert payload["metrics"] is None
    assert payload["output_path"] == str(output_csv)
    extras = payload.get("extras")
    assert isinstance(extras, dict)
    assert "rows" in extras and "processed" in extras

    assert {path for path, _ in written} == {output_csv}
    assert captured.activities == [("ACT1", "ACT2", "ACT3")]


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_activity_pipeline__happy_path(
    activity_resource_dir: Path, cfg, tmp_path, monkeypatch
):
    _configure_cfg(cfg)
    input_csv = _copy_resource(activity_resource_dir, "ids_happy.csv", tmp_path)
    output_csv = tmp_path / "activities.csv"
    chunk_df = pd.read_csv(activity_resource_dir / "chunk_happy.csv")

    monkeypatch.setattr(
        get_activity_data,
        "_load_assay_src_lookup",
        lambda *_: {
            "ASSAY1": "SRC-ASSAY1",
            "ASSAY2": "SRC-ASSAY2",
            "ASSAY3": "SRC-ASSAY3",
        },
    )

    captured = _install_fetch_stubs(monkeypatch, chunk_df)
    written = _install_writer_stub(monkeypatch)
    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr("library.validation.logger", logger_stub)

    args = _make_args(input_csv, output_csv)

    exit_code = get_activity_data.run_chembl(cfg, args)

    activity_events = [event for _, event, _ in logger_stub.events]

    assert exit_code == 0
    assert {path for path, _ in written} == {output_csv}
    assert set(activity_events).issuperset(
        {
            "activity_pipeline_start",
            "activity_http_config",
            "activity_pipeline_done",
            "records_dropped",
            "schema_validate_start",
            "schema_validate_done",
        }
    )
    config_payloads = [
        payload
        for level, event, payload in logger_stub.events
        if event == "activity_http_config"
    ]
    assert config_payloads
    config_payload = config_payloads[0]
    assert config_payload["activity_batch_size"] == cfg.activity.batch_size
    assert config_payload["activity_timeout"] == cfg.activity.timeout
    assert config_payload["api_retries"] == cfg.api.retries
    assert config_payload["retry_max_attempts"] == cfg.retry.max_attempts
    assert captured.activities == [("ACT1", "ACT2", "ACT3")]
    assert captured.testitems

    written_df = written[0][1]
    assert list(written_df["activity_id"]) == ["ACT1", "ACT2", "ACT3"]
    assert pd.api.types.is_float_dtype(written_df["standard_value"])  # type: ignore[arg-type]
    assert written_df["standard_value"].tolist() == [5.5, 7.25, 9.0]
    assert "src_assay_id" in written_df.columns
    src_assay_series = written_df["src_assay_id"].astype("string")
    assert src_assay_series.tolist() == ["SRC-ASSAY1", "SRC-ASSAY2", "SRC-ASSAY3"]
    assert src_assay_series.str.strip().ne("").all()
    completion_payloads = [
        payload
        for _, event, payload in logger_stub.events
        if event == "activity_pipeline_completion"
    ]
    assert completion_payloads
    summary_payload = completion_payloads[-1]
    assert summary_payload["mode"] == "run"
    assert summary_payload["rows"] == 3
    total_cells = written_df.size
    expected_null_fraction = (
        float(np.count_nonzero(written_df.isna().to_numpy()) / total_cells)
        if total_cells
        else 0.0
    )
    assert summary_payload["null_fraction"] == pytest.approx(expected_null_fraction)
    assert summary_payload["output"] == str(output_csv)
    report_path = cfg.io.output_dir / "activity.postprocess.report.json"
    assert report_path.exists()
    report_payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert report_payload["table"] == "activity"
    assert "metrics" in report_payload

    meta_path = output_csv.with_name(output_csv.name + ".meta.yaml")
    assert meta_path.exists()
    meta = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    dictionaries = meta.get("dictionaries")
    assert isinstance(dictionaries, dict)

    root_resource = get_resource("dictionary_root")
    target_resource = get_resource("target_types")

    assert dictionaries.get("dictionary_root") == {
        "version": root_resource.version,
        "sha256": root_resource.sha256,
    }
    assert dictionaries.get("target_types") == {
        "version": target_resource.version,
        "sha256": target_resource.sha256,
    }


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_activity_pipeline__extended_columns_dtype_coercion(
    activity_resource_dir: Path, cfg, tmp_path, monkeypatch
):
    _configure_cfg(cfg)
    input_csv = _copy_resource(activity_resource_dir, "ids_happy.csv", tmp_path)
    output_csv = tmp_path / "activities.csv"

    chunk_df = pd.read_csv(activity_resource_dir / "chunk_happy.csv")
    chunk_df["activity_chembl_id"] = pd.Series(
        [float("nan"), 502.0, float("nan")], dtype="float64"
    )
    chunk_df["compound_name"] = pd.Series(
        [float("nan"), 1.0, float("nan")], dtype="float64"
    )
    chunk_df["molecule_pref_name"] = pd.Series(
        ["NALTREXONE", "CLOZAPINE", pd.NA], dtype="string"
    )
    chunk_df["pchembl_value"] = pd.Series(
        ["7.5", "not-a-number", "6.25"], dtype="object"
    )
    chunk_df["log_value"] = pd.Series([float("nan")] * len(chunk_df), dtype="float64")

    testitem_frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
            "pref_name": ["NALTREXONE", "CLOZAPINE", "OXFORD"],
        }
    )

    monkeypatch.setattr(
        get_activity_data,
        "_load_assay_src_lookup",
        lambda *_: {
            "ASSAY1": "SRC-ASSAY1",
            "ASSAY2": "SRC-ASSAY2",
            "ASSAY3": "SRC-ASSAY3",
        },
    )

    _install_fetch_stubs(monkeypatch, chunk_df, testitem_frame=testitem_frame)
    written = _install_writer_stub(monkeypatch)
    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr("library.validation.logger", logger_stub)

    args = _make_args(input_csv, output_csv)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        exit_code = get_activity_data.run_chembl(cfg, args)

    dtype_warnings = [
        warning
        for warning in caught
        if "incompatible dtype" in str(warning.message)
        or "Setting an item of incompatible dtype" in str(warning.message)
    ]

    assert exit_code == 0
    assert not dtype_warnings
    assert {path for path, _ in written} == {output_csv}

    debug_events = [event for level, event, _ in logger_stub.events if level == "debug"]
    assert "activity_extended_fallback_conversion_failed" in debug_events

    written_df = written[0][1]
    assert pd.api.types.is_string_dtype(written_df["activity_chembl_id"])  # type: ignore[arg-type]
    assert pd.api.types.is_string_dtype(written_df["compound_name"])  # type: ignore[arg-type]
    assert pd.api.types.is_float_dtype(written_df["log_value"])  # type: ignore[arg-type]
    assert written_df.loc[0, "activity_chembl_id"] == "ACT1"
    assert written_df.loc[0, "compound_name"] == "NALTREXONE"
    assert pytest.approx(written_df.loc[0, "log_value"], rel=1e-6) == 7.5
    assert pd.isna(written_df.loc[1, "log_value"])
    assert written_df["standard_value"].tolist() == [5.5, 7.25, 9.0]
    assert "src_assay_id" in written_df.columns
    src_assay_series = written_df["src_assay_id"].astype("string")
    assert src_assay_series.tolist() == ["SRC-ASSAY1", "SRC-ASSAY2", "SRC-ASSAY3"]
    assert src_assay_series.str.strip().ne("").all()


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_activity_pipeline__missing_column_input(
    activity_resource_dir: Path, cfg, tmp_path, monkeypatch
):
    _configure_cfg(cfg)
    input_csv = _copy_resource(
        activity_resource_dir, "ids_missing_column.csv", tmp_path
    )
    output_csv = tmp_path / "activities.csv"

    # Ensure we never reach the fetch stage when validation of inputs fails.
    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        _DummyChemblClient,
    )
    monkeypatch.setattr(
        get_activity_data.cl, "get_activities", lambda *_, **__: pd.DataFrame()
    )
    monkeypatch.setattr(get_activity_data, "_load_assay_src_lookup", lambda *_: {})
    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr("library.validation.logger", logger_stub)

    args = _make_args(input_csv, output_csv)

    exit_code = get_activity_data.run_chembl(cfg, args)

    activity_events = [event for _, event, _ in logger_stub.events]

    assert exit_code == 1
    assert "read_fail" in activity_events
    assert not output_csv.exists()


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
@pytest.mark.parametrize("runner_variant", ["cli", "api"])
def test_activity_pipeline__batch_size_clamped(
    cfg,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    runner_variant,
) -> None:
    _configure_cfg(cfg)
    cfg.activity.batch_size = get_activity_data.MAX_ACTIVITY_CHUNK_SIZE + 5

    input_csv = tmp_path / "ids.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    output_csv = tmp_path / "activities.csv"

    logger_stub = _RecordingLogger()
    _patch_activity_loggers(monkeypatch, logger_stub)

    captured_batch_sizes: list[int | None] = []

    def _fake_run_chembl(config: object, args: argparse.Namespace) -> int:
        captured_batch_sizes.append(getattr(config.activity, "batch_size", None))
        return 0

    args = _make_args(input_csv, output_csv)

    exit_code = _invoke_activity_runner(
        cfg,
        args,
        runner_variant,
        monkeypatch,
        runner=_fake_run_chembl,
    )

    assert exit_code == 0
    assert captured_batch_sizes == [get_activity_data.MAX_ACTIVITY_CHUNK_SIZE]

    warning_events = [
        (level, event, payload)
        for level, event, payload in logger_stub.events
        if level == "warning"
    ]
    assert any(event == "activity_batch_size_clamped" for _, event, _ in warning_events)


def test_activity_pipeline__malformed_values(
    activity_resource_dir: Path, cfg, tmp_path, monkeypatch
):
    _configure_cfg(cfg)
    input_csv = _copy_resource(activity_resource_dir, "ids_happy.csv", tmp_path)
    output_csv = tmp_path / "activities.csv"
    chunk_df = pd.read_csv(activity_resource_dir / "chunk_malformed.csv")

    monkeypatch.setattr(get_activity_data, "_load_assay_src_lookup", lambda *_: {})
    _install_fetch_stubs(monkeypatch, chunk_df)
    written = _install_writer_stub(monkeypatch)
    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr("library.validation.logger", logger_stub)

    args = _make_args(input_csv, output_csv)

    exit_code = get_activity_data.run_chembl(cfg, args)

    activity_events = [event for _, event, _ in logger_stub.events]

    assert exit_code == 1
    assert not written
    assert "validation_failed" in activity_events
    failure_path = output_csv.with_name("activities_failure_cases.csv")
    assert failure_path.exists()
    failure_df = pd.read_csv(failure_path)
    assert len(failure_df) >= 2
    assert set(failure_df["column"]) == {"standard_value"}
    failure_cases = " ".join(str(value) for value in failure_df["failure_case"])
    assert "not-a-number" in failure_cases
    assert ">=" in failure_cases
    error_events = {event for level, event, _ in logger_stub.events if level == "error"}
    assert {"validation_failed", "activity_pipeline_failed"}.issubset(error_events)


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_activity_pipeline__deduplicates_identifiers(
    activity_resource_dir: Path, cfg, tmp_path, monkeypatch
):
    _configure_cfg(cfg)
    input_csv = _copy_resource(activity_resource_dir, "ids_duplicates.csv", tmp_path)
    output_csv = tmp_path / "activities.csv"
    chunk_df = pd.read_csv(activity_resource_dir / "chunk_happy.csv")

    monkeypatch.setattr(get_activity_data, "_load_assay_src_lookup", lambda *_: {})
    captured = _install_fetch_stubs(monkeypatch, chunk_df)
    written = _install_writer_stub(monkeypatch)
    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr("library.validation.logger", logger_stub)

    args = _make_args(input_csv, output_csv)

    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 0
    assert len(captured.activities) == 1
    recorded = captured.activities[0]
    assert len({item.lower() for item in recorded}) == 3
    assert len(written) == 1
    written_df = written[0][1]
    assert written_df["activity_id"].tolist() == ["ACT1", "ACT2", "ACT3"]
    info_events = {event for _, event, _ in logger_stub.events}
    assert "records_dropped" in info_events
    assert "schema_validate_done" in info_events
    warning_events = {
        event for level, event, _ in logger_stub.events if level == "warning"
    }
    assert "read_ids_dropped_na_markers" not in warning_events


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_activity_pipeline__fills_compound_name_from_pref_name(
    cfg, tmp_path, monkeypatch
):
    _configure_cfg(cfg)

    total_records = 40
    filled_records = 38  # 95%

    input_csv = tmp_path / "ids_many.csv"
    records = [f"ACT{i}" for i in range(1, total_records + 1)]
    input_csv.write_text("activity_id\n" + "\n".join(records) + "\n", encoding="utf-8")

    chunk_records = []
    for idx, activity in enumerate(records, start=1):
        chunk_records.append(
            {
                "activity_id": activity,
                "molecule_chembl_id": f"CHEMBL{idx}",
                "assay_chembl_id": f"ASSAY{idx}",
                "standard_value": float(idx),
                "standard_units": "nM",
                "standard_type": "IC50",
                "relation": "=",
            }
        )
    chunk_df = pd.DataFrame.from_records(chunk_records)

    pref_name_records = []
    for idx in range(1, filled_records + 1):
        pref_name_records.append(
            {
                "molecule_chembl_id": f"CHEMBL{idx}",
                "pref_name": f"Compound {idx}",
            }
        )
    testitem_df = pd.DataFrame.from_records(pref_name_records)

    monkeypatch.setattr(get_activity_data, "_load_assay_src_lookup", lambda *_: {})
    capture = _install_fetch_stubs(monkeypatch, chunk_df, testitem_frame=testitem_df)
    written = _install_writer_stub(monkeypatch)
    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr("library.validation.logger", logger_stub)

    output_csv = tmp_path / "activities.csv"
    args = _make_args(input_csv, output_csv)

    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 0
    assert capture.testitems, "expected test item enrichment to be invoked"
    assert written, "pipeline should write output"

    written_df = written[0][1]
    compound_series = written_df["compound_name"].astype("string")
    fill_mask = compound_series.notna() & compound_series.str.strip().ne("")
    fill_rate = fill_mask.sum() / len(compound_series)
    assert fill_rate >= 0.95


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
@pytest.mark.parametrize(
    ("error_factory", "expected_message"),
    [
        pytest.param(
            lambda: requests.Timeout("pref name lookup timed out"),
            "pref name lookup timed out",
            id="timeout",
        ),
        pytest.param(
            lambda: requests.HTTPError("404 Client Error: Not Found for url"),
            "404 Client Error",
            id="http-404",
        ),
    ],
)
def test_activity_pipeline__records_pref_name_fetch_failures(
    cfg,
    tmp_path,
    monkeypatch,
    error_factory,
    expected_message,
):
    _configure_cfg(cfg)

    input_csv = tmp_path / "ids_pref_name.csv"
    input_csv.write_text("activity_id\nACT1\nACT2\n", encoding="utf-8")

    chunk_df = pd.DataFrame.from_records(
        [
            {
                "activity_id": "ACT1",
                "molecule_chembl_id": "CHEMBL1",
                "assay_chembl_id": "ASSAY1",
                "standard_value": 1.0,
                "standard_units": "nM",
                "standard_type": "IC50",
                "relation": "=",
            },
            {
                "activity_id": "ACT2",
                "molecule_chembl_id": "CHEMBL2",
                "assay_chembl_id": "ASSAY2",
                "standard_value": 2.0,
                "standard_units": "nM",
                "standard_type": "IC50",
                "relation": "=",
            },
        ]
    )

    capture = _install_fetch_stubs(
        monkeypatch,
        chunk_df,
        testitem_error=error_factory,
    )
    written = _install_writer_stub(monkeypatch)
    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)
    monkeypatch.setattr("library.validation.logger", logger_stub)

    output_csv = tmp_path / "activities.csv"
    args = _make_args(input_csv, output_csv)

    exit_code = get_activity_data.run_chembl(cfg, args)

    assert exit_code == 0
    assert capture.testitems, "expected pref name fetch attempts"
    assert written, "pipeline should produce output"

    pref_name_events = [
        payload
        for level, event, payload in logger_stub.events
        if level == "warning" and event == "pref_name_fetch_failed"
    ]
    assert pref_name_events, "pref name fetch failure should be logged"
    assert pref_name_events[0]["pending"] == list(capture.testitems[0])
    assert expected_message in str(pref_name_events[0]["error"])

    fetch_failure_path = output_csv.with_name("activities_fetch_failures.csv")
    assert fetch_failure_path.exists()

    failure_df = pd.read_csv(fetch_failure_path)
    assert len(failure_df) == 1
    recorded_ids = failure_df.loc[0, "chunk_ids"].split(",")
    assert recorded_ids == list(capture.testitems[0])
    assert failure_df.loc[0, "chunk_size"] == len(capture.testitems[0])
    assert expected_message in str(failure_df.loc[0, "error"])

    meta_path = Path(f"{fetch_failure_path}.meta.yaml")
    assert meta_path.exists()
    meta = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert set(meta.get("columns", [])) >= {"chunk_ids", "chunk_size", "error"}


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_load_assay_src_lookup__real_dictionary_includes_known_assays() -> None:
    lookup = get_activity_data._load_assay_src_lookup(DICTIONARY_DIR)

    assert lookup.get("CHEMBL1762864") == "357280"
    assert lookup.get("CHEMBL1762866") == "357284"
    assert all(str(key).strip() and str(value).strip() for key, value in lookup.items())
