from __future__ import annotations

import argparse
import io
from pathlib import Path
from typing import Sequence

import pytest

from scripts import get_data
from tests.helpers.logs import iter_events, parse_log_lines


def _make_config(tmp_path: Path) -> get_data.PipelineRunConfig:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    for name, filename in get_data._DEFAULT_INPUT_FILES.items():
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
    )


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
    steps = get_data._PIPELINE_STEPS
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
def test_run_pipeline__propagates_step_failure(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = _make_config(tmp_path)
    failure_calls: list[Sequence[str]] = []

    def _success(argv: Sequence[str]) -> int:
        final_out = Path(argv[argv.index("--final-out") + 1])
        final_out.parent.mkdir(parents=True, exist_ok=True)
        final_out.write_text("id\n1\n", encoding="utf-8")
        return 0

    def _failure(argv: Sequence[str]) -> int:
        failure_calls.append(tuple(argv))
        return 2

    steps = (
        get_data.PipelineStep("document", _success, None),
        get_data.PipelineStep("target", _failure, None),
        get_data.PipelineStep("assay", _success, None),
    )

    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="unit")
    )
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", steps, raising=False)
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    status = get_data.run_pipeline(cfg)
    assert status == 2
    assert failure_calls, "expected failing step to execute"
    records = parse_log_lines(stream.getvalue())
    events = list(iter_events(records))
    assert "step_failed" in events
    assert events[-1] != "step_done"


@pytest.mark.unit
def test_run_pipeline__handles_step_exception(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = _make_config(tmp_path)

    def _raising(argv: Sequence[str]) -> int:  # pragma: no cover - executed via pipeline
        raise RuntimeError("malformed output")

    steps = (get_data.PipelineStep("document", _raising, None),)
    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="unit")
    )
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", steps, raising=False)
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    status = get_data.run_pipeline(cfg)
    assert status == 1
    records = parse_log_lines(stream.getvalue())
    assert any(entry["event"] == "step_exception" for entry in records)
