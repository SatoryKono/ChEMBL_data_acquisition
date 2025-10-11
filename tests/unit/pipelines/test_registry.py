from __future__ import annotations

from pathlib import Path

import pytest

from library.pipelines.registry import PipelineStep, load_pipeline_registry


@pytest.mark.unit
def test_default_registry__exposes_expected_steps() -> None:
    steps = load_pipeline_registry()
    names = [step.name for step in steps]
    assert names == ["document", "target", "assay", "testitem", "activity"]

    document = _step_by_name(steps, "document")
    assert document.extra_args == ("--mode", "all")
    assert document.input_filename == "document.csv"
    assert document.output_stem == "documents"
    assert document.produces == ("documents",)
    assert document.consumes == ()

    target = _step_by_name(steps, "target")
    assert target.subcommand == "all"
    assert target.input_filename == "target.csv"
    assert target.output_stem == "targets"
    assert target.produces == ("targets",)
    assert target.consumes == ()

    assay = _step_by_name(steps, "assay")
    assert assay.subcommand is None
    assert assay.supports_dry_run is False
    assert assay.produces == ("assays",)

    testitem = _step_by_name(steps, "testitem")
    assert callable(testitem.main)
    assert testitem.produces == ("testitems",)

    activity = _step_by_name(steps, "activity")
    assert activity.supports_dry_run is True
    assert activity.input_filename == "activity.csv"
    assert activity.output_stem == "activity"
    assert activity.produces == ("activity", "activities")
    assert activity.consumes == ("documents", "targets", "assays", "testitems")


@pytest.mark.unit
def test_registry_loader__supports_yaml_overrides(tmp_path: Path) -> None:
    registry_path = tmp_path / "pipelines.yaml"
    registry_path.write_text(
        """
        pipelines:
          - name: example
            callable: scripts.get_document_data:main
            input: example.csv
            output: examples
            flags:
              dry_run: true
        """.strip()
        + "\n",
        encoding="utf-8",
    )

    steps = load_pipeline_registry(registry_path)
    assert [step.name for step in steps] == ["example"]
    [step] = steps
    assert step.input_filename == "example.csv"
    assert step.output_stem == "examples"
    assert step.supports_dry_run is True
    assert step.produces == ("examples",)


def _step_by_name(steps: tuple[PipelineStep, ...], name: str) -> PipelineStep:
    for step in steps:
        if step.name == name:
            return step
    raise AssertionError(f"step {name!r} not found")
