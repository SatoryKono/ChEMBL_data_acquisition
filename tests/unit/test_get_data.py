from __future__ import annotations

import argparse
import io
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from typing import Sequence

import pytest

from scripts import get_data
from tests.helpers.logs import iter_events, parse_log_lines
from tests.helpers.manifests import load_latest_manifest
from library.pipelines.common import PipelineRunResult


def _make_config(tmp_path: Path) -> get_data.PipelineRunConfig:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    input_files = dict(get_data.DEFAULT_INPUT_FILES)
    output_stems = dict(get_data.DEFAULT_OUTPUT_STEMS)
    for name, filename in input_files.items():
        target = input_dir / filename
        target.write_text("id\nplaceholder\n", encoding="utf-8")
    config_path = base_path / "config.yaml"
    config_path.write_text("io:\n  csv_sep: ','\n", encoding="utf-8")
    return get_data.PipelineRunConfig(
        base_path=base_path,
        input_dir=input_dir,
        output_dir=output_dir,
        config_path=config_path,
        date_prefix="20200101",
        log_level="INFO",
        limit=None,
        force=False,
        skip_existing=False,
        dry_run=False,
        input_files=input_files,
        output_stems=output_stems,
    )


def _load_manifest(cfg: get_data.PipelineRunConfig) -> dict[str, object]:
    _, manifest = load_latest_manifest(cfg.base_path)
    return manifest


@pytest.mark.unit
def test_parse_args__defaults(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("CHEMBL_DA_BASE_PATH", raising=False)
    args = get_data._parse_args([])
    assert args.base_path == Path("data")
    assert args.input_dir == Path("input")
    assert args.output_dir == Path("output")
    assert args.config == get_data.DEFAULT_CONFIG_PATH
    expected_prefix = get_data.datetime.now(get_data.UTC).strftime("%Y%m%d")
    assert args.date_prefix == expected_prefix
    assert args.log_level == "INFO"
    assert args.verbose is False
    assert args.limit is None
    assert args.force is False
    assert args.skip_existing is False
    assert args.dry_run is False
    assert args.pipeline_registry is None
    assert args.override_input == []
    assert args.override_output_stem == []
    assert args.override_subcommand == []


@pytest.mark.unit
def test_parse_args__custom_paths(tmp_path: Path) -> None:
    base = tmp_path / "workspace"
    args = get_data._parse_args(
        [
            "--base-path",
            str(base),
            "--input-dir",
            "inputs",
            "--output-dir",
            "outputs",
            "--config",
            str(tmp_path / "custom.yaml"),
            "--date",
            "20240203",
            "--log-level",
            "debug",
            "--limit",
            "50",
            "--force",
            "--skip-existing",
            "--dry-run",
            "--verbose",
            "--pipeline-registry",
            str(tmp_path / "registry.yaml"),
            "--override-input",
            "document=document_custom.csv",
            "--override-output-stem",
            "target=custom_targets",
            "--override-subcommand",
            "target=sync",
        ]
    )
    assert args.base_path == base
    assert args.input_dir == Path("inputs")
    assert args.output_dir == Path("outputs")
    assert args.config == tmp_path / "custom.yaml"
    assert args.date_prefix == "20240203"
    assert args.log_level == "debug"
    assert args.verbose is True
    assert args.limit == 50
    assert args.force is True
    assert args.skip_existing is True
    assert args.dry_run is True
    assert args.pipeline_registry == tmp_path / "registry.yaml"
    assert args.override_input == ["document=document_custom.csv"]
    assert args.override_output_stem == ["target=custom_targets"]
    assert args.override_subcommand == ["target=sync"]


@pytest.mark.unit
def test_prepare_config__verbose_overrides_level(tmp_path: Path) -> None:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    config_path = base_path / "config.yaml"
    config_path.write_text("io:\n  csv_sep: ','\n", encoding="utf-8")

    args = argparse.Namespace(
        base_path=base_path,
        input_dir=Path("input"),
        output_dir=Path("output"),
        config=config_path,
        date_prefix="20240204",
        log_level="info",
        limit=None,
        force=False,
        skip_existing=False,
        dry_run=False,
        verbose=True,
    )

    cfg = get_data._prepare_config(args)
    assert cfg.log_level == "DEBUG"


@pytest.mark.unit
def test_pipeline_step_registration__expected_shape() -> None:
    steps = get_data.DEFAULT_PIPELINE_STEPS
    names = [step.name for step in steps]
    assert names == ["document", "target", "assay", "testitem", "activity"]
    assert steps[0].extra_args == ("--mode", "all")
    assert steps[1].subcommand == "all"
    assert steps[-1].supports_dry_run is True
    for step in steps:
        assert callable(step.main)


@pytest.mark.unit
def test_configure_logging__delegates_to_configure_logger(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: list[get_data.LoggerConfig] = []
    original_configure = get_data.configure_logger

    def _wrapper(cfg: get_data.LoggerConfig) -> get_data.Logger:
        captured.append(cfg)
        return original_configure(cfg)

    monkeypatch.setattr(get_data, "configure_logger", _wrapper)
    logger = get_data._configure_logging("warn", run_id="fixed")
    assert captured, "expected configure_logger to be invoked"
    cfg = captured[0]
    assert cfg.level == "WARN"
    assert cfg.run_id == "fixed"


@pytest.mark.unit
def test_configure_logging__invalid_level() -> None:
    with pytest.raises(ValueError):
        get_data._configure_logging("invalid")


@pytest.mark.unit
def test_resolve_pipeline_steps__applies_overrides() -> None:
    args = argparse.Namespace(
        pipeline_registry=None,
        override_input=["document=doc.csv"],
        override_output_stem=["activity=custom"],
        override_subcommand=["target="],
    )

    steps = get_data._resolve_pipeline_steps(args)
    mapping = {step.name: step for step in steps}
    assert mapping["document"].input_filename == "doc.csv"
    assert mapping["activity"].output_stem == "custom"
    assert mapping["target"].subcommand is None


@pytest.mark.unit
def test_run_pipeline__propagates_step_failure(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = _make_config(tmp_path)
    failure_calls: list[Sequence[str]] = []

    def _success_runner(_: get_data.Config, options: SimpleNamespace) -> PipelineRunResult:
        path = Path(options.output_csv)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("id\n1\n", encoding="utf-8")
        return PipelineRunResult(exit_code=0, output_path=path, executed=True, written=True)

    def _failure_runner(_: get_data.Config, options: SimpleNamespace) -> PipelineRunResult:
        failure_calls.append((str(options.input_csv), str(options.output_csv)))
        return PipelineRunResult(exit_code=2, output_path=Path(options.output_csv), executed=True)

    def _build_options(
        cfg: get_data.PipelineRunConfig, input_path: Path, output_path: Path
    ) -> SimpleNamespace:
        return SimpleNamespace(input_csv=input_path, output_csv=output_path)

    steps = (
        get_data.PipelineStep(
            name="document",
            main=lambda _: 0,
            input_filename="document.csv",
            output_stem="documents",
        ),
        get_data.PipelineStep(
            name="target",
            main=lambda _: 0,
            input_filename="target.csv",
            output_stem="targets",
        ),
        get_data.PipelineStep(
            name="assay",
            main=lambda _: 0,
            input_filename="assay.csv",
            output_stem="assays",
        ),
    )

    apis = {
        "document": get_data.PipelineApi(_build_options, _success_runner),
        "target": get_data.PipelineApi(_build_options, _failure_runner),
        "assay": get_data.PipelineApi(_build_options, _success_runner),
    }

    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="unit")
    )
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    monkeypatch.setattr(get_data, "_PIPELINE_APIS", apis, raising=False)
    status = get_data.run_pipeline(cfg, steps=steps)
    assert status == 2
    assert failure_calls, "expected failing step to execute"

    manifest = _load_manifest(cfg)
    assert manifest["run"]["exit_code"] == 2
    step_entries = manifest["steps"]
    assert len(step_entries) == 3
    first, second, third = step_entries
    assert first["status"] == "success"
    assert first["output"]["exists"] is True
    assert first["output"]["checksum_sha256"]
    assert second["status"] == "failed"
    assert second["executed"] is True
    assert second["exit_code"] == 2
    assert second["reason"] == "non_zero_exit"
    assert second["output"]["exists"] is False
    assert third["status"] == "pending"

    records = parse_log_lines(stream.getvalue())
    events = list(iter_events(records))
    assert "step_failed" in events
    assert events[-1] != "step_done"


@pytest.mark.unit
def test_run_pipeline__handles_step_exception(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = _make_config(tmp_path)

    def _raising_runner(_: get_data.Config, options: SimpleNamespace) -> PipelineRunResult:
        raise RuntimeError("malformed output")

    def _build_options(
        cfg: get_data.PipelineRunConfig, input_path: Path, output_path: Path
    ) -> SimpleNamespace:
        return SimpleNamespace(input_csv=input_path, output_csv=output_path)

    steps = (
        get_data.PipelineStep(
            name="document",
            main=lambda _: 0,
            input_filename="document.csv",
            output_stem="documents",
        ),
    )
    apis = {"document": get_data.PipelineApi(_build_options, _raising_runner)}
    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="unit")
    )
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    monkeypatch.setattr(get_data, "_PIPELINE_APIS", apis, raising=False)
    status = get_data.run_pipeline(cfg, steps=steps)
    assert status == 1
    manifest = _load_manifest(cfg)
    assert manifest["run"]["exit_code"] == 1
    step_entries = manifest["steps"]
    assert len(step_entries) == 1
    step_entry = step_entries[0]
    assert step_entry["status"] == "failed"
    assert step_entry["reason"] == "exception"
    assert step_entry["output"]["exists"] is False
    records = parse_log_lines(stream.getvalue())
    assert any(entry["event"] == "step_exception" for entry in records)


@pytest.mark.unit
def test_run_pipeline__dry_run_manifest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = _make_config(tmp_path)
    cfg = replace(cfg, dry_run=True)

    executions: list[int] = []

    def _record(argv: Sequence[str]) -> int:
        executions.append(len(argv))
        return 0

    steps = (
        get_data.PipelineStep(
            "document",
            _record,
            "document.csv",
            "documents",
        ),
        get_data.PipelineStep(
            "target",
            _record,
            "target.csv",
            "targets",
        ),
    )

    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="unit"),
    )
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", steps, raising=False)
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    status = get_data.run_pipeline(cfg, steps=steps)
    assert status == 0
    assert not executions

    manifest = _load_manifest(cfg)
    assert manifest["run"]["dry_run"] is True
    assert manifest["run"]["exit_code"] == 0
    step_entries = manifest["steps"]
    assert [entry["status"] for entry in step_entries] == ["skipped", "skipped"]
    for entry in step_entries:
        assert entry["reason"] == "dry_run"
        assert entry["output"]["exists"] is False
