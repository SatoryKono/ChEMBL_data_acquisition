"""Integration tests for the activity pipeline CLI helpers."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path
from typing import Iterable

import pandas as pd
import pytest

from dataclasses import dataclass

from scripts import get_activity_data


class _DummyChemblClient:
    """Minimal stand-in for :class:`ChemblClient` used in tests."""

    def __init__(self, *args, **kwargs) -> None:  # pragma: no cover - interface compatibility
        pass

    def __enter__(self) -> "_DummyChemblClient":  # pragma: no cover - trivial
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:  # pragma: no cover - trivial
        return False


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
    target.write_text((resource_dir / name).read_text(encoding="utf-8"), encoding="utf-8")
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
        if identifiers:
            mask = testitem_frame["molecule_chembl_id"].astype(str).isin(identifiers)
            return testitem_frame.loc[mask].copy().reset_index(drop=True)
        return testitem_frame.iloc[0:0].copy()

    monkeypatch.setattr(get_activity_data.cl, "get_testitem", _fake_get_testitem)
    monkeypatch.setattr(get_activity_data, "ChemblClient", _DummyChemblClient)
    return _FetchCapture(captured_activities, captured_testitems)


def _install_writer_stub(monkeypatch: pytest.MonkeyPatch) -> list[tuple[Path, pd.DataFrame]]:
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
            order = [str(col) for col in col_order]
            result = result.reindex(columns=order, fill_value=pd.NA)
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


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_activity_pipeline__timeout_clamped_when_below_minimum(cfg, tmp_path, monkeypatch):
    _configure_cfg(cfg)
    cfg.activity.timeout = get_activity_data.MIN_ACTIVITY_TIMEOUT - 5

    input_csv = tmp_path / "ids.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    output_csv = tmp_path / "activities.csv"

    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)

    captured_timeout: dict[str, float] = {}

    def _run_chembl_stub(passed_cfg, _passed_args):
        captured_timeout["timeout"] = float(passed_cfg.activity.timeout)
        return 0

    monkeypatch.setattr(get_activity_data, "run_chembl", _run_chembl_stub)

    args = _make_args(input_csv, output_csv)

    exit_code = get_activity_data.run(cfg, args)

    warning_events = [event for level, event, _ in logger_stub.events if level == "warning"]
    assert "activity_timeout_clamped" in warning_events
    assert captured_timeout["timeout"] == pytest.approx(get_activity_data.MIN_ACTIVITY_TIMEOUT)
    assert cfg.activity.timeout == pytest.approx(get_activity_data.MIN_ACTIVITY_TIMEOUT)
    assert exit_code == 0


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_activity_pipeline__warns_when_retry_disabled(cfg, tmp_path, monkeypatch):
    _configure_cfg(cfg)
    cfg.retry.max_attempts = 1

    input_csv = tmp_path / "ids.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    output_csv = tmp_path / "activities.csv"

    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)

    def _run_stub(passed_cfg, _args):
        assert passed_cfg.retry.max_attempts == 1
        return 0

    monkeypatch.setattr(get_activity_data, "run_chembl", _run_stub)

    args = _make_args(input_csv, output_csv)

    exit_code = get_activity_data.run(cfg, args)

    warning_events = [event for level, event, _ in logger_stub.events if level == "warning"]
    assert "activity_retry_disabled" in warning_events
    assert exit_code == 0


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_activity_pipeline__warns_when_api_retries_disabled(cfg, tmp_path, monkeypatch):
    _configure_cfg(cfg)
    cfg.api.retries = 0

    input_csv = tmp_path / "ids.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    output_csv = tmp_path / "activities.csv"

    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)

    def _run_stub(passed_cfg, _args):
        assert passed_cfg.api.retries == 0
        return 0

    monkeypatch.setattr(get_activity_data, "run_chembl", _run_stub)

    args = _make_args(input_csv, output_csv)

    exit_code = get_activity_data.run(cfg, args)

    warning_events = [event for level, event, _ in logger_stub.events if level == "warning"]
    assert "activity_api_retry_disabled" in warning_events
    assert exit_code == 0
@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_activity_pipeline__happy_path(activity_resource_dir: Path, cfg, tmp_path, monkeypatch):
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
        payload for level, event, payload in logger_stub.events if event == "activity_http_config"
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
    completion_messages = [
        event
        for _, event, _ in logger_stub.events
        if event.startswith("Completed get_activity_data pipeline:")
    ]
    assert completion_messages
    summary_message = completion_messages[-1]
    assert "mode=run" in summary_message
    assert "rows=3" in summary_message


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_activity_pipeline__extended_columns_dtype_coercion(
    activity_resource_dir: Path, cfg, tmp_path, monkeypatch
):
    _configure_cfg(cfg)
    input_csv = _copy_resource(activity_resource_dir, "ids_happy.csv", tmp_path)
    output_csv = tmp_path / "activities.csv"

    chunk_df = pd.read_csv(activity_resource_dir / "chunk_happy.csv")
    chunk_df["activity_chembl_id"] = pd.Series([float("nan"), 502.0, float("nan")], dtype="float64")
    chunk_df["compound_name"] = pd.Series([float("nan"), 1.0, float("nan")], dtype="float64")
    chunk_df["molecule_pref_name"] = pd.Series(["NALTREXONE", "CLOZAPINE", pd.NA], dtype="string")
    chunk_df["pchembl_value"] = pd.Series(["7.5", "not-a-number", "6.25"], dtype="object")
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
        lambda *_: {"ASSAY1": "SRC-ASSAY1", "ASSAY2": "SRC-ASSAY2", "ASSAY3": "SRC-ASSAY3"},
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
def test_activity_pipeline__extended_columns_dtype_coercion(
    activity_resource_dir: Path, cfg, tmp_path, monkeypatch
):
    _configure_cfg(cfg)
    input_csv = _copy_resource(activity_resource_dir, "ids_happy.csv", tmp_path)
    output_csv = tmp_path / "activities.csv"

    chunk_df = pd.read_csv(activity_resource_dir / "chunk_happy.csv")
    chunk_df["activity_chembl_id"] = pd.Series([float("nan"), 502.0, float("nan")], dtype="float64")
    chunk_df["compound_name"] = pd.Series([float("nan"), 1.0, float("nan")], dtype="float64")
    chunk_df["molecule_pref_name"] = pd.Series(["NALTREXONE", "CLOZAPINE", pd.NA], dtype="string")
    chunk_df["pchembl_value"] = pd.Series(["7.5", "not-a-number", "6.25"], dtype="object")
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
        lambda *_: {"ASSAY1": "SRC-ASSAY1", "ASSAY2": "SRC-ASSAY2", "ASSAY3": "SRC-ASSAY3"},
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
def test_activity_pipeline__missing_column_input(activity_resource_dir: Path, cfg, tmp_path, monkeypatch):
    _configure_cfg(cfg)
    input_csv = _copy_resource(activity_resource_dir, "ids_missing_column.csv", tmp_path)
    output_csv = tmp_path / "activities.csv"

    # Ensure we never reach the fetch stage when validation of inputs fails.
    monkeypatch.setattr(get_activity_data, "ChemblClient", _DummyChemblClient)
    monkeypatch.setattr(get_activity_data.cl, "get_activities", lambda *_, **__: pd.DataFrame())
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
def test_activity_pipeline__batch_size_clamped(cfg, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _configure_cfg(cfg)
    cfg.activity.batch_size = get_activity_data.MAX_ACTIVITY_CHUNK_SIZE + 5

    input_csv = tmp_path / "ids.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    output_csv = tmp_path / "activities.csv"

    logger_stub = _RecordingLogger()
    monkeypatch.setattr(get_activity_data, "logger", logger_stub)

    captured_batch_sizes: list[int | None] = []

    def _fake_run_chembl(config: object, args: argparse.Namespace) -> int:
        captured_batch_sizes.append(getattr(config.activity, "batch_size", None))
        return 0

    monkeypatch.setattr(get_activity_data, "run_chembl", _fake_run_chembl)

    args = _make_args(input_csv, output_csv)

    exit_code = get_activity_data.run(cfg, args)

    assert exit_code == 0
    assert captured_batch_sizes == [get_activity_data.MAX_ACTIVITY_CHUNK_SIZE]

    warning_events = [
        (level, event, payload)
        for level, event, payload in logger_stub.events
        if level == "warning"
    ]
    assert any(event == "activity_batch_size_clamped" for _, event, _ in warning_events)

def test_activity_pipeline__malformed_values(activity_resource_dir: Path, cfg, tmp_path, monkeypatch):
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
def test_activity_pipeline__deduplicates_identifiers(activity_resource_dir: Path, cfg, tmp_path, monkeypatch):
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
    warning_events = {event for level, event, _ in logger_stub.events if level == "warning"}
    assert "read_ids_dropped_na_markers" not in warning_events


@pytest.mark.integration
@pytest.mark.usefixtures("deterministic_env")
def test_activity_pipeline__fills_compound_name_from_pref_name(cfg, tmp_path, monkeypatch):
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
