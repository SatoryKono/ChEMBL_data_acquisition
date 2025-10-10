from __future__ import annotations

import argparse
import io
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from library.pipelines.common import PipelineRunResult
from scripts import get_data
from tests.helpers.logs import iter_events, parse_log_lines
from tests.helpers.manifests import load_latest_manifest


def _make_config(tmp_path: Path) -> get_data.PipelineRunConfig:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    input_files = dict(get_data.DEFAULT_INPUT_FILES)
    output_stems = dict(get_data.DEFAULT_OUTPUT_STEMS)
    subcommands = {
        step.name: step.subcommand for step in get_data.DEFAULT_PIPELINE_STEPS
    }
    for _name, filename in input_files.items():
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
        subcommands=subcommands,
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
    assert args.print_config is False
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
            "--print-config",
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
    assert args.print_config is True
    assert args.pipeline_registry == tmp_path / "registry.yaml"
    assert args.override_input == ["document=document_custom.csv"]
    assert args.override_output_stem == ["target=custom_targets"]
    assert args.override_subcommand == ["target=sync"]


@pytest.mark.unit
def test_main__print_config__exits_early(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base_path = tmp_path
    input_dir = base_path / "input"
    input_dir.mkdir()
    (base_path / "output").mkdir()
    config_path = base_path / "config.yaml"
    config_path.write_text("io:\n  csv_sep: ','\n", encoding="utf-8")

    printed: dict[str, object] = {}

    def _fake_print_config(cfg: object) -> None:
        printed["config"] = cfg

    monkeypatch.setattr(get_data, "print_config", _fake_print_config)

    loaded: dict[str, Path] = {}

    def _fake_load_config(path: Path, *, base_path: Path) -> object:
        loaded["path"] = path
        loaded["base"] = base_path
        return object()

    monkeypatch.setattr(get_data, "load_config", _fake_load_config)

    run_invoked = False

    def _fake_run_pipeline(
        cfg: get_data.PipelineRunConfig,
        *,
        steps: Sequence[get_data.PipelineStep] | None = None,
    ) -> int:
        nonlocal run_invoked
        run_invoked = True
        return 0

    monkeypatch.setattr(get_data, "run_pipeline", _fake_run_pipeline)
    monkeypatch.setattr(get_data, "_resolve_pipeline_steps", lambda args=None: ())

    def _fake_setup_cli_logging(
        script_name: str,
        log_cfg: get_data.LoggerConfig,
        date_str: str | None = None,
        *,
        log_dir: Path | None = None,
    ) -> object:
        stream = io.StringIO()
        ctx_cfg = get_data.LoggerConfig(
            level=log_cfg.level,
            run_id=log_cfg.run_id,
            stream=stream,
            handlers=list(log_cfg.handlers),
            redact_secrets=log_cfg.redact_secrets,
            logger_name=log_cfg.logger_name,
        )

        class _Manager:
            def __enter__(self) -> SimpleNamespace:
                return SimpleNamespace(log_cfg=ctx_cfg, console_stream=stream)

            def __exit__(self, exc_type, exc, tb) -> bool:
                return False

        return _Manager()

    monkeypatch.setattr(get_data, "setup_cli_logging", _fake_setup_cli_logging)

    status = get_data.main(
        [
            "--base-path",
            str(base_path),
            "--input-dir",
            "input",
            "--output-dir",
            "output",
            "--config",
            str(config_path),
            "--print-config",
        ]
    )

    assert status == 0
    assert loaded["path"] == config_path.resolve()
    assert loaded["base"] == base_path.resolve()
    assert printed["config"] is not None
    assert run_invoked is False


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
def test_configure_logging__delegates_to_configure_logger(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: list[get_data.LoggerConfig] = []
    original_configure = get_data.configure_logger

    def _wrapper(cfg: get_data.LoggerConfig) -> get_data.Logger:
        captured.append(cfg)
        return original_configure(cfg)

    monkeypatch.setattr(get_data, "configure_logger", _wrapper)
    get_data._configure_logging("warn", run_id="fixed")
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
def test_override_subcommand__target_pipeline_uses_selected_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    config_path = base_path / "config.yaml"
    config_path.write_text("{}", encoding="utf-8")
    (input_dir / "target.csv").write_text("target_chembl_id\nT1\n", encoding="utf-8")

    args = argparse.Namespace(
        base_path=base_path,
        input_dir=Path("input"),
        output_dir=Path("output"),
        config=config_path,
        date_prefix="20240214",
        log_level="INFO",
        limit=None,
        force=False,
        skip_existing=False,
        dry_run=False,
        verbose=False,
        pipeline_registry=None,
        override_input=[],
        override_output_stem=[],
        override_subcommand=["target=chembl"],
    )

    resolved_steps = get_data._resolve_pipeline_steps(args)
    target_only = tuple(step for step in resolved_steps if step.name == "target")
    cfg = get_data._prepare_config(args, target_only)
    assert cfg.subcommand_for("target") == "chembl"

    captured: list[str] = []

    def _build_options(
        run_cfg: get_data.PipelineRunConfig, input_path: Path, output_path: Path
    ) -> SimpleNamespace:
        return SimpleNamespace(
            input_csv=input_path,
            output_csv=output_path,
            command=run_cfg.subcommand_for("target") or "all",
        )

    def _runner(_: get_data.Config, options: SimpleNamespace) -> PipelineRunResult:
        captured.append(options.command)
        destination = Path(options.output_csv)
        destination.write_text("target_chembl_id\nT1\n", encoding="utf-8")
        return PipelineRunResult(
            exit_code=0,
            output_path=destination,
            executed=True,
            reason=None,
            written=True,
        )

    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="unit_override")
    )
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)
    monkeypatch.setattr(
        get_data,
        "_PIPELINE_APIS",
        {"target": get_data.PipelineApi(_build_options, _runner)},
        raising=False,
    )
    monkeypatch.setattr(
        get_data,
        "load_config",
        lambda *args, **kwargs: SimpleNamespace(),
        raising=False,
    )

    status = get_data.run_pipeline(cfg, steps=target_only)
    assert status == 0
    assert captured == ["chembl"]


@pytest.mark.unit
def test_run_pipeline__propagates_step_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = _make_config(tmp_path)
    failure_calls: list[Sequence[str]] = []

    def _success_runner(
        _: get_data.Config, options: SimpleNamespace
    ) -> PipelineRunResult:
        path = Path(options.output_csv)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("id\n1\n", encoding="utf-8")
        return PipelineRunResult(
            exit_code=0, output_path=path, executed=True, written=True
        )

    def _failure_runner(
        _: get_data.Config, options: SimpleNamespace
    ) -> PipelineRunResult:
        failure_calls.append((str(options.input_csv), str(options.output_csv)))
        return PipelineRunResult(
            exit_code=2, output_path=Path(options.output_csv), executed=True
        )

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
def test_run_pipeline__handles_step_exception(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = _make_config(tmp_path)

    def _raising_runner(
        _: get_data.Config, options: SimpleNamespace
    ) -> PipelineRunResult:
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
def test_run_pipeline__dry_run_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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
