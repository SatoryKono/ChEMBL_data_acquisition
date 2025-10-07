from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from config import CONFIG_SCHEMA_PATH
from library.config.postprocessing import (
    PostprocessingConfigError,
    diff_postprocessing_pipeline,
    list_postprocessing_tables,
    load_postprocessing_pipeline,
    main,
)


def _write(path: Path, content: str) -> None:
    path.write_text(textwrap.dedent(content).strip() + "\n", encoding="utf-8")


def _prepare_config(tmp_path: Path) -> Path:
    config_dir = tmp_path / "postprocessing"
    config_dir.mkdir()

    _write(
        config_dir / "defaults.yaml",
        """
        version: 1
        defaults:
          step:
            enabled: true
            continue_on_error: false
            kwargs:
              columns: []
        flags:
          qa_checks: true
        steps:
          - id: ensure-columns
            callable: library.postprocessing.pipeline:passthrough
            kwargs:
              columns: []
          - id: run-qa
            callable: library.postprocessing.pipeline:passthrough
            gates:
              - qa_checks
        """,
    )

    _write(
        config_dir / "target.yaml",
        """
        version: 2
        flags:
          qa_checks: false
          include_cellularity: true
        steps:
          - id: ensure-columns
            kwargs:
              columns:
                - pref_name
          - id: finalize
            description: "Finalize export"
            callable: library.postprocessing.pipeline:passthrough
            gates:
              - include_cellularity
        """,
    )
    return config_dir


@pytest.mark.unit
def test_load_postprocessing_pipeline__merges_defaults_and_overrides(tmp_path: Path) -> None:
    config_dir = _prepare_config(tmp_path)

    pipeline = load_postprocessing_pipeline(
        "target", base_dir=config_dir, schema_path=CONFIG_SCHEMA_PATH
    )

    assert pipeline.version == 2
    assert pipeline.flags == {"qa_checks": False, "include_cellularity": True}
    assert [step.identifier for step in pipeline.steps] == [
        "ensure-columns",
        "run-qa",
        "finalize",
    ]
    ensure_step, qa_step, finalize_step = pipeline.steps
    assert ensure_step.kwargs["columns"] == ["pref_name"]
    assert qa_step.gates == ("qa_checks",)
    assert finalize_step.gates == ("include_cellularity",)


@pytest.mark.unit
def test_load_postprocessing_pipeline__raises_for_unknown_gate(tmp_path: Path) -> None:
    config_dir = _prepare_config(tmp_path)
    _write(
        config_dir / "target.yaml",
        """
        steps:
          - id: ensure-columns
            callable: library.postprocessing.pipeline:passthrough
            gates:
              - undefined_gate
        """,
    )

    with pytest.raises(PostprocessingConfigError) as excinfo:
        load_postprocessing_pipeline(
            "target", base_dir=config_dir, schema_path=CONFIG_SCHEMA_PATH
        )

    assert "unknown gate" in str(excinfo.value)


@pytest.mark.unit
def test_load_postprocessing_pipeline__rejects_new_step_without_callable(tmp_path: Path) -> None:
    config_dir = _prepare_config(tmp_path)
    _write(
        config_dir / "target.yaml",
        """
        steps:
          - id: ensure-columns
            kwargs:
              columns:
                - pref_name
          - id: missing-callable
        """,
    )

    with pytest.raises(PostprocessingConfigError) as excinfo:
        load_postprocessing_pipeline(
            "target", base_dir=config_dir, schema_path=CONFIG_SCHEMA_PATH
        )

    assert "must define a callable" in str(excinfo.value)


@pytest.mark.unit
def test_postprocessing_config_cli__validate_and_diff(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    config_dir = _prepare_config(tmp_path)

    exit_code = main(
        [
            "--config-dir",
            str(config_dir),
            "--schema",
            str(CONFIG_SCHEMA_PATH),
            "validate",
        ]
    )
    captured = capsys.readouterr()
    assert exit_code == 0
    assert "Validated post-processing configs: target" in captured.out

    exit_code = main(
        [
            "--config-dir",
            str(config_dir),
            "--schema",
            str(CONFIG_SCHEMA_PATH),
            "diff",
            "target",
        ]
    )
    diff_output = capsys.readouterr().out
    assert exit_code == 0
    assert "--- defaults" in diff_output
    assert "+++ target" in diff_output


@pytest.mark.unit
def test_list_postprocessing_tables__discovers_yaml_files(tmp_path: Path) -> None:
    config_dir = _prepare_config(tmp_path)
    (config_dir / "extra.txt").write_text("ignored", encoding="utf-8")

    tables = list_postprocessing_tables(config_dir)

    assert tables == ["target"]


@pytest.mark.unit
def test_diff_postprocessing_pipeline__empty_when_matching(tmp_path: Path) -> None:
    config_dir = _prepare_config(tmp_path)
    _write(
        config_dir / "target.yaml",
        """
        version: 1
        flags:
          qa_checks: true
        steps:
          - id: ensure-columns
            callable: library.postprocessing.pipeline:passthrough
            kwargs:
              columns: []
          - id: run-qa
            callable: library.postprocessing.pipeline:passthrough
            gates:
              - qa_checks
        """,
    )

    diff_output = diff_postprocessing_pipeline(
        "target", base_dir=config_dir, schema_path=CONFIG_SCHEMA_PATH
    )

    assert diff_output.strip() == ""
