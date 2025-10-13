from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from library.pipelines.common import PipelineRunResult
from library.project_version import get_pipeline_version
from scripts import get_data
from tests.helpers.manifests import list_manifest_files, load_latest_manifest

_TEST_GIT_SHA = "scheduler-test-sha"


def _build_run_config(
    tmp_path: Path, steps: tuple[get_data.PipelineStep, ...]
) -> get_data.PipelineRunConfig:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    config_path = base_path / "config.yaml"
    config_path.write_text("{}", encoding="utf-8")
    subcommands = get_data.PipelineSubcommands.from_mapping(
        {step.name: step.subcommand for step in steps}
    )
    return get_data.PipelineRunConfig(
        base_path=base_path,
        input_dir=input_dir,
        output_dir=output_dir,
        config_path=config_path,
        date_prefix="20240101",
        log_level="INFO",
        limit=None,
        force=False,
        skip_existing=False,
        dry_run=False,
        rerun_postprocess=False,
        input_files=get_data.PipelineInputFiles.from_mapping(
            {step.name: step.input_filename for step in steps}
        ),
        output_stems=get_data.PipelineOutputStems.from_mapping(
            {step.name: step.output_stem for step in steps}
        ),
        subcommands=subcommands,
    )


def _make_pipeline_api(name: str, executed: list[str]) -> get_data.PipelineApi:
    def _builder(
        cfg: get_data.PipelineRunConfig, input_path: Path, output_path: Path
    ) -> SimpleNamespace:
        return SimpleNamespace(input_csv=input_path, output_csv=output_path)

    def _runner(config: object, options: SimpleNamespace) -> PipelineRunResult:
        executed.append(name)
        output_path = Path(options.output_csv)
        output_path.write_text("value\n1\n", encoding="utf-8")
        return PipelineRunResult(
            exit_code=0,
            output_path=output_path,
            executed=True,
            reason=None,
            written=True,
        )

    return get_data.PipelineApi(_builder, _runner)


@pytest.mark.e2e
@pytest.mark.smoke
def test_scheduler__selective_run_respects_dependencies(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    steps = (
        get_data.PipelineStep(
            name="audit",
            main=lambda _: 0,
            input_filename="audit.csv",
            output_stem="audit",
            depends_on=("transform",),
            produces=("audit",),
            consumes=("transformed",),
        ),
        get_data.PipelineStep(
            name="transform",
            main=lambda _: 0,
            input_filename="transform.csv",
            output_stem="transformed",
            produces=("transformed",),
            consumes=("prep",),
        ),
    )

    cfg = _build_run_config(tmp_path, steps)
    (cfg.input_dir / "transform.csv").write_text("id\n1\n", encoding="utf-8")
    (cfg.input_dir / "audit.csv").write_text("id\n1\n", encoding="utf-8")
    external = cfg.output_dir / f"output.prep_{cfg.date_prefix}.csv"
    external.write_text("id\n1\n", encoding="utf-8")

    executed: list[str] = []
    monkeypatch.setattr(
        get_data,
        "_PIPELINE_APIS",
        {step.name: _make_pipeline_api(step.name, executed) for step in steps},
        raising=False,
    )
    monkeypatch.setattr(
        get_data,
        "load_config",
        lambda *args, **kwargs: SimpleNamespace(),
        raising=False,
    )
    monkeypatch.setenv("GIT_SHA", _TEST_GIT_SHA)

    status = get_data.run_pipeline(cfg, steps=steps)
    assert status == 0
    assert executed == ["transform", "audit"]

    transform_output = cfg.output_dir / f"output.transformed_{cfg.date_prefix}.csv"
    audit_output = cfg.output_dir / f"output.audit_{cfg.date_prefix}.csv"
    assert transform_output.exists()
    assert audit_output.exists()

    manifests = list_manifest_files(cfg.base_path)
    assert len(manifests) == 1
    _, manifest = load_latest_manifest(cfg.base_path)
    step_names = [entry["name"] for entry in manifest["steps"]]
    assert step_names == ["transform", "audit"]
    assert manifest["run"]["status"] == "success"
    run_info = manifest["run"]
    assert run_info["pipeline_version"] == get_pipeline_version()
    assert run_info["git_sha"] == _TEST_GIT_SHA


@pytest.mark.e2e
def test_scheduler__fails_on_missing_external_dependency(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    steps = (
        get_data.PipelineStep(
            name="audit",
            main=lambda _: 0,
            input_filename="audit.csv",
            output_stem="audit",
            depends_on=("transform",),
            produces=("audit",),
            consumes=("transformed",),
        ),
        get_data.PipelineStep(
            name="transform",
            main=lambda _: 0,
            input_filename="transform.csv",
            output_stem="transformed",
            produces=("transformed",),
            consumes=("prep",),
        ),
    )

    cfg = _build_run_config(tmp_path, steps)
    (cfg.input_dir / "transform.csv").write_text("id\n1\n", encoding="utf-8")
    (cfg.input_dir / "audit.csv").write_text("id\n1\n", encoding="utf-8")

    executed: list[str] = []
    monkeypatch.setattr(
        get_data,
        "_PIPELINE_APIS",
        {step.name: _make_pipeline_api(step.name, executed) for step in steps},
        raising=False,
    )
    monkeypatch.setattr(
        get_data,
        "load_config",
        lambda *args, **kwargs: SimpleNamespace(),
        raising=False,
    )
    monkeypatch.setenv("GIT_SHA", _TEST_GIT_SHA)

    status = get_data.run_pipeline(cfg, steps=steps)
    assert status == 1
    assert executed == []

    manifests = list_manifest_files(cfg.base_path)
    assert len(manifests) == 1
    _, manifest = load_latest_manifest(cfg.base_path)
    assert manifest["run"]["status"] == "failed"
    run_info = manifest["run"]
    assert run_info["pipeline_version"] == get_pipeline_version()
    assert run_info["git_sha"] == _TEST_GIT_SHA
    transform_entry = next(
        entry for entry in manifest["steps"] if entry["name"] == "transform"
    )
    assert transform_entry["status"] == "failed"
    assert transform_entry["reason"] == "dependency_missing"
    assert transform_entry["executed"] is False
