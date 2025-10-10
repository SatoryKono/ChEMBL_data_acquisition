"""Unit tests for :mod:`scripts.get_testitem_data`."""

from __future__ import annotations

import argparse
from collections.abc import Iterable, MutableMapping
from pathlib import Path

import pandas as pd
import pytest

from library.config import Config
from scripts import get_testitem_data


class _MemoryLogger:
    """Capture structured log events emitted by the pipeline module."""

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
    monkeypatch.setattr(get_testitem_data, "logger", logger)
    return logger


@pytest.mark.parametrize(
    "values, expected",
    [
        (["chembl1", " CHEMBL2 "], ["CHEMBL1", "CHEMBL2"]),
        ([None, ""], ["", ""]),
        (["ChEmBl3"], ["CHEMBL3"]),
        (["chembl4", " chembl5\t"], ["CHEMBL4", "CHEMBL5"]),
        (["Σ"], ["Σ"]),
        (["  "], [""]),
    ],
)
def test_normalise_chembl_ids__cases(
    values: Iterable[str | None], expected: list[str]
) -> None:
    series = pd.Series(values, dtype="object")

    result = get_testitem_data._normalise_chembl_ids(series)

    assert result.tolist() == expected
    assert str(result.dtype) == "string"


def test_load_molecule_hierarchy_lookup__delegates_and_logs(
    cfg: Config,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping = {"chembl1": "chembl_parent", "chembl2": None}
    normalised: list[tuple[str, str | None]] = []

    def fake_loader(path: str, encoding: str, delimiter: str) -> dict[str, str | None]:
        assert path.endswith("hierarchy.csv")
        assert encoding == cfg.io.csv_encoding
        assert delimiter == cfg.io.csv_sep
        return mapping

    def fake_normalise(value: str | None, *, child_id: str) -> str | None:
        normalised.append((child_id, value))
        return value.upper() if value else None

    monkeypatch.setattr(
        get_testitem_data, "_load_molecule_hierarchy_mapping", fake_loader
    )
    monkeypatch.setattr(
        get_testitem_data, "_normalise_parent_identifier", fake_normalise
    )

    lookup = get_testitem_data.load_molecule_hierarchy_lookup(
        Path("hierarchy.csv"), io_cfg=cfg.io
    )

    assert lookup == {"chembl1": "CHEMBL_PARENT", "chembl2": None}
    assert normalised == [
        ("chembl1", "chembl_parent"),
        ("chembl2", None),
    ]
    assert (
        "info",
        "molecule_hierarchy_lookup_loaded",
        {"path": "hierarchy.csv", "rows": 1},
    ) in logger_stub.events


def test_load_molecule_hierarchy_lookup__missing_file_returns_empty(
    cfg: Config,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_loader(path: str, encoding: str, delimiter: str) -> dict[str, str | None]:
        raise FileNotFoundError(path)

    monkeypatch.setattr(
        get_testitem_data, "_load_molecule_hierarchy_mapping", fake_loader
    )

    lookup = get_testitem_data.load_molecule_hierarchy_lookup(
        Path("missing.csv"), io_cfg=cfg.io
    )

    assert lookup == {}
    assert (
        "info",
        "molecule_hierarchy_lookup_missing",
        {"path": "missing.csv"},
    ) in logger_stub.events


def test_load_molecule_hierarchy_lookup__invalid_data_raises(
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_loader(path: str, encoding: str, delimiter: str) -> dict[str, str | None]:
        raise ValueError("bad delimiter")

    monkeypatch.setattr(
        get_testitem_data, "_load_molecule_hierarchy_mapping", fake_loader
    )

    with pytest.raises(ValueError, match="invalid hierarchy lookup: bad delimiter"):
        get_testitem_data.load_molecule_hierarchy_lookup(Path("bad.csv"), io_cfg=cfg.io)


def test_run_chembl__passes_options(
    cfg: Config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    captured: dict[str, object] = {}

    def fake_run_pipeline(
        config: Config, options: get_testitem_data.TestitemPipelineOptions
    ) -> int:
        captured["config"] = config
        captured["options"] = options
        return 0

    monkeypatch.setattr(get_testitem_data, "run_testitem_pipeline", fake_run_pipeline)

    args = argparse.Namespace(input_csv=input_csv, final_out=output_csv)

    exit_code = get_testitem_data.run_chembl(cfg, args)

    assert exit_code == 0
    options = captured["options"]
    assert isinstance(options, get_testitem_data.TestitemPipelineOptions)
    assert options.input_csv == input_csv
    assert options.output_csv == output_csv


def test_run_chembl__passes_limit_and_offset(
    cfg: Config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\nCHEMBL2\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    captured: dict[str, object] = {}

    def fake_run_pipeline(
        config: Config, options: get_testitem_data.TestitemPipelineOptions
    ) -> int:
        captured["options"] = options
        return 0

    monkeypatch.setattr(get_testitem_data, "run_testitem_pipeline", fake_run_pipeline)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        limit=5,
        offset=2,
    )

    exit_code = get_testitem_data.run_chembl(cfg, args)

    assert exit_code == 0
    options = captured["options"]
    assert isinstance(options, get_testitem_data.TestitemPipelineOptions)
    assert options.limit == 5
    assert options.offset == 2


def test_run__skip_existing_without_force(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"
    output_csv.write_text("existing", encoding="utf-8")

    called = False

    def fake_run(cfg: Config, args: argparse.Namespace) -> int:
        nonlocal called
        called = True
        return 0

    monkeypatch.setattr(get_testitem_data, "run_chembl", fake_run)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=True,
        force=False,
    )

    exit_code = get_testitem_data.run(cfg, args)

    assert exit_code == 0
    assert not called
    assert (
        "info",
        "pipeline_skip_existing",
        {"output": str(output_csv)},
    ) in logger_stub.events


def test_run__force_overrides_skip(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"
    output_csv.write_text("existing", encoding="utf-8")

    calls: list[str] = []

    def fake_run(cfg: Config, args: argparse.Namespace) -> int:
        calls.append("run_chembl")
        return 0

    monkeypatch.setattr(get_testitem_data, "run_chembl", fake_run)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=True,
        force=True,
    )

    exit_code = get_testitem_data.run(cfg, args)

    assert exit_code == 0
    assert calls == ["run_chembl"]


def test_run__non_zero_exit_logs_error(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    monkeypatch.setattr(get_testitem_data, "run_chembl", lambda *_: 5)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
    )

    exit_code = get_testitem_data.run(cfg, args)

    assert exit_code == 5
    assert (
        "error",
        "testitem_pipeline_failed",
        {"output": str(output_csv), "exit_code": 5},
    ) in logger_stub.events


def test_run__success_logs_completion(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    monkeypatch.setattr(get_testitem_data, "run_chembl", lambda *_: 0)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
    )

    exit_code = get_testitem_data.run(cfg, args)

    assert exit_code == 0
    assert (
        "info",
        "Postprocessing skipped (flag --postprocess not set)",
        {},
    ) in logger_stub.events
    assert (
        "info",
        "testitem_pipeline_done",
        {"output": str(output_csv)},
    ) in logger_stub.events


def test_run__postprocess_enabled_runs_pipeline(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.testitems.csv"
    output_csv.write_text("id,name\nCHEMBL1,Example\n", encoding="utf-8")

    postprocessed_csv = output_csv.with_name("output_postprocessed.testitems.csv")
    report_path = tmp_path / "testitems.postprocess.report.json"

    class DummyMetrics:
        pipeline_version = "v1"
        validation = type("Validation", (), {"schema": "TestitemPostprocess"})()

        def summary(self) -> dict[str, object]:
            return {
                "pipeline_version": self.pipeline_version,
                "rows": 1,
                "columns": 2,
                "duration_s": 0.5,
                "steps": 3,
            }

    dummy_report_path = report_path

    class DummyResult:
        metrics = DummyMetrics()
        output_path = postprocessed_csv
        report_path = dummy_report_path

    pipeline_cfg = type("PipelineCfg", (), {"pipeline_version": "v1"})()
    csv_cfg = type("CsvCfg", (), {})()
    captured: dict[str, object] = {}

    def fake_get_pipeline_config(table: str, override: object) -> object:
        captured["pipeline_table"] = table
        captured["pipeline_override"] = override
        return pipeline_cfg

    def fake_get_csv_runtime_config(config: object) -> object:
        assert config is pipeline_cfg
        return csv_cfg

    def fake_run_postprocessing_pipeline(
        table: str,
        input_path: Path,
        destination: Path,
        runtime_cfg: object,
    ) -> DummyResult:
        captured.update(
            postprocess_table=table,
            postprocess_input=input_path,
            postprocess_output=destination,
            runtime_cfg=runtime_cfg,
        )
        return DummyResult()

    def fake_runner(frame, *, pipeline_version: str | None = None, logger=None):  # type: ignore[no-untyped-def]
        return frame, DummyMetrics()

    def fake_validator(frame):  # type: ignore[no-untyped-def]
        return frame

    import library.postprocess as postprocess_mod
    import library.postprocessing.testitem as postprocess_testitem

    monkeypatch.setattr(postprocess_mod, "get_pipeline_config", fake_get_pipeline_config)
    monkeypatch.setattr(postprocess_mod, "get_csv_runtime_config", fake_get_csv_runtime_config)
    monkeypatch.setattr(
        postprocess_mod,
        "run_postprocessing_pipeline",
        fake_run_postprocessing_pipeline,
    )
    monkeypatch.setattr(postprocess_testitem, "run_testitem_pipeline", fake_runner)
    monkeypatch.setattr(postprocess_testitem, "validate_testitems", fake_validator)
    monkeypatch.setattr(postprocess_testitem, "TESTITEM_SCHEMA", object())

    monkeypatch.setattr(get_testitem_data, "run_chembl", lambda *_: 0)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
        postprocess=True,
        config=None,
    )

    exit_code = get_testitem_data.run(cfg, args)

    assert exit_code == 0
    expected_destination = postprocessed_csv
    assert captured["postprocess_table"] == "testitems"
    assert captured["postprocess_input"] == output_csv
    assert captured["postprocess_output"] == expected_destination
    assert captured["runtime_cfg"].pipeline_config is pipeline_cfg
    assert captured["runtime_cfg"].csv_runtime_config is csv_cfg

    assert (
        "info",
        "testitem_pipeline_done",
        {
            "output": str(output_csv),
            "postprocess_output": str(postprocessed_csv),
            "postprocess_pipeline_version": "v1",
            "postprocess_rows": 1,
            "postprocess_columns": 2,
            "postprocess_duration_s": 0.5,
            "postprocess_steps": 3,
            "postprocess_schema": "TestitemPostprocess",
            "postprocess_report": str(report_path),
        },
    ) in logger_stub.events


def test_add_pubchem_data__delegates_to_pipeline(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_add_pubchem(
        df: pd.DataFrame,
        cfg: get_testitem_data.PubChemCfg,
        *,
        client: object | None,
        api_cfg: object | None,
        timeout: float | None,
        cid_cache: MutableMapping[str, str | None] | None,
        resolution_cache: MutableMapping[tuple[str | None, ...], object] | None,
        parent_record_cache: MutableMapping[str, pd.Series | None] | None,
        testitem_fields: Iterable[str] | None,
        request_limit: int,
    ) -> pd.DataFrame:
        captured.update(
            df=df,
            cfg=cfg,
            client=client,
            api_cfg=api_cfg,
            timeout=timeout,
            cid_cache=cid_cache,
            resolution_cache=resolution_cache,
            parent_record_cache=parent_record_cache,
            testitem_fields=tuple(testitem_fields or ()),
            request_limit=request_limit,
        )
        return df.assign(attached=True)

    monkeypatch.setattr(
        get_testitem_data.pipeline, "add_pubchem_data", fake_add_pubchem
    )

    frame = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    cfg = get_testitem_data.PubChemCfg()
    cid_cache: dict[str, str | None] = {}
    resolution_cache: dict[tuple[str | None, ...], object] = {}
    parent_cache: dict[str, pd.Series | None] = {}

    result = get_testitem_data.add_pubchem_data(
        frame,
        cfg,
        client=None,
        api_cfg=None,
        timeout=3.0,
        cid_cache=cid_cache,
        resolution_cache=resolution_cache,
        parent_record_cache=parent_cache,
        testitem_fields=["molecule_chembl_id"],
        request_limit=10,
    )

    assert result["attached"].tolist() == [True]
    assert captured["cid_cache"] is cid_cache
    assert captured["resolution_cache"] is resolution_cache
    assert captured["parent_record_cache"] is parent_cache
    assert captured["testitem_fields"] == ("molecule_chembl_id",)
