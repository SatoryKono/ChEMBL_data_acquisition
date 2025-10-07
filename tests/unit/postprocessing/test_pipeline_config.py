from __future__ import annotations

from pathlib import Path

import textwrap

import pytest

from library.postprocess.config import PipelineConfigError, PipelineStep, load_pipeline_config


def test_load_pipeline_config__expands_base_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    base_dir = tmp_path / "data"
    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(base_dir))

    cfg = load_pipeline_config("activities")

    assert cfg.pipeline_version
    assert any(step.name == "compute_bounds" for step in cfg.steps)

    compute_step = next(step for step in cfg.steps if step.name == "compute_bounds")
    input_path = compute_step.params["input_csv"]
    expected = (base_dir / "output" / "activities" / "output.activity.csv").resolve()
    assert Path(input_path).resolve() == expected


def test_pipeline_step__resolves_callable(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(tmp_path))

    cfg = load_pipeline_config("assays")
    step = cfg.steps[0]
    callable_obj = step.resolve()
    result = callable_obj(**dict(step.params))

    assert isinstance(result, dict)
    assert result["step"].startswith("assays.")
    assert set(result["params"]).issuperset(step.params.keys())


def test_load_pipeline_config__skips_disabled_steps(tmp_path: Path) -> None:
    config_text = textwrap.dedent(
        """
        pipeline_version: "test"
        steps:
          - name: active
            enabled: true
            callable: "library.postprocess.activities.steps:compute_bounds"
            params:
              input_csv: "in.csv"
              rounding_digits: 2
              clamp_nonnegative: false
          - name: disabled
            enabled: false
            callable: "library.postprocess.activities.steps:annotate_targets"
            params:
              dictionary_csv: "targets.csv"
              fallback_label: "unknown"
              minimum_confidence: 0.5
        """
    ).strip()
    config_path = tmp_path / "pipeline.yaml"
    config_path.write_text(config_text, encoding="utf-8")

    cfg = load_pipeline_config("activities", path=config_path)

    assert [step.name for step in cfg.steps] == ["active"]

    with pytest.raises(PipelineConfigError):
        PipelineStep(name="broken", callable_path="not:callable", params={}).resolve()
