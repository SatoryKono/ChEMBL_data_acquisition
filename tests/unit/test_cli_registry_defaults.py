"""Tests ensuring CLI defaults track the pipeline registry."""

from __future__ import annotations

from pathlib import Path

import pytest

import library.pipelines.registry as registry_module
from library.pipelines.registry import PipelineStep
from scripts import get_assay_data


def _patch_registry(monkeypatch: pytest.MonkeyPatch, step: PipelineStep) -> None:
    def _fake_load_pipeline_registry(
        source: object | None = None,
    ) -> tuple[PipelineStep, ...]:
        return (step,)

    monkeypatch.setattr(
        registry_module, "load_pipeline_registry", _fake_load_pipeline_registry
    )


def _make_assay_step(
    *,
    input_filename: str = "registry_assay.csv",
    output_stem: str = "registry_assays",
) -> PipelineStep:
    return PipelineStep(
        name="assay",
        main=get_assay_data.main,
        input_filename=input_filename,
        output_stem=output_stem,
        subcommand=None,
        output_flag="--final-out",
        extra_args=(),
        supports_dry_run=False,
        callable_path="scripts.get_assay_data:main",
    )


@pytest.mark.unit
def test_assay_cli_defaults__parser_uses_registry(monkeypatch: pytest.MonkeyPatch) -> None:
    step = _make_assay_step(input_filename="patched_assay.csv")
    _patch_registry(monkeypatch, step)

    parser, _ = get_assay_data.build_parser()
    args = parser.parse_args([])

    assert args.input_csv == Path(step.input_filename)


@pytest.mark.unit
def test_assay_cli_defaults__main_uses_registry(monkeypatch: pytest.MonkeyPatch) -> None:
    step = _make_assay_step(
        input_filename="patched_assay_cli.csv",
        output_stem="patched_assay_output",
    )
    _patch_registry(monkeypatch, step)

    captured: dict[str, object] = {}
    original_prepare = get_assay_data.cli.prepare_io_paths

    def _prepare_io_paths(
        args,
        *,
        input_default: str | Path | None = None,
        output_stem: str | None = None,
        **kwargs,
    ) -> None:
        captured["input_default"] = input_default
        captured["output_stem"] = output_stem
        return original_prepare(
            args,
            input_default=input_default,
            output_stem=output_stem,
            **kwargs,
        )

    monkeypatch.setattr(get_assay_data.cli, "prepare_io_paths", _prepare_io_paths)
    monkeypatch.setattr(get_assay_data, "run_cli_command", lambda **_: 0)

    class _DummyLoggingContext:
        def __init__(self) -> None:
            self.log_cfg = get_assay_data.LoggerConfig()

        def __enter__(self) -> "_DummyLoggingContext":
            return self

        def __exit__(self, exc_type, exc, tb) -> bool:
            return False

    monkeypatch.setattr(
        get_assay_data, "setup_cli_logging", lambda *args, **kwargs: _DummyLoggingContext()
    )

    exit_code = get_assay_data.main([])

    assert exit_code == 0
    assert captured["input_default"] == step.input_filename
    assert captured["output_stem"] == step.output_stem
