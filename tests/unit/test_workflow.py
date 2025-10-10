from dataclasses import dataclass
from pathlib import Path

import pytest

from library.orchestration.workflow import (
    PreparedPipelineStep,
    StepExecutionResult,
    execute_workflow,
    temporary_output_path,
)
from library.pipelines.registry import PipelineStep


@dataclass
class _ConfigStub:
    base_path: Path
    input_dir: Path
    output_dir: Path
    config_path: Path
    date_prefix: str
    log_level: str
    limit: int | None
    force: bool
    skip_existing: bool
    dry_run: bool
    input_files: dict[str, str]
    output_stems: dict[str, str]

    def input_path(self, name: str) -> Path:
        return self.input_dir / self.input_files[name]

    def output_path(self, name: str) -> Path:
        stem = self.output_stems[name]
        return self.output_dir / f"output.{stem}_{self.date_prefix}.csv"


def _make_config(tmp_path: Path) -> _ConfigStub:
    base_path = tmp_path
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    config_path = tmp_path / "config.yaml"
    input_dir.mkdir()
    output_dir.mkdir()
    config_path.write_text("{}", encoding="utf-8")
    return _ConfigStub(
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
        input_files={"doc": "document.csv", "target": "target.csv"},
        output_stems={"doc": "documents", "target": "targets"},
    )


@pytest.mark.unit
def test_execute_workflow__passes_resolved_paths(tmp_path: Path) -> None:
    cfg = _make_config(tmp_path)
    (cfg.input_dir / "document.csv").write_text("id\n1\n", encoding="utf-8")

    step = PipelineStep(
        name="doc",
        main=lambda args: 0,
        input_filename="document.csv",
        output_stem="documents",
    )

    calls: list[tuple[_ConfigStub, Path, Path, Path]] = []

    def runner(
        cfg_arg: _ConfigStub,
        input_path: Path,
        final_output: Path,
        working_output: Path,
    ) -> StepExecutionResult:
        calls.append((cfg_arg, input_path, final_output, working_output))
        return StepExecutionResult(exit_code=0, executed=True, status="success")

    prepared = PreparedPipelineStep(step=step, invoke=runner)

    results = list(execute_workflow(cfg, [prepared]))

    assert len(results) == 1
    [call] = calls
    assert call[0] is cfg
    assert call[1] == cfg.input_dir / "document.csv"
    expected_output = cfg.output_dir / "output.documents_20240101.csv"
    assert call[2] == expected_output
    assert call[3] == temporary_output_path(expected_output)
    result = results[0]
    assert result.step is step
    assert result.result.exit_code == 0
    assert result.result.status == "success"


@pytest.mark.unit
def test_execute_workflow__stops_after_failure(tmp_path: Path) -> None:
    cfg = _make_config(tmp_path)
    (cfg.input_dir / "document.csv").write_text("id\n1\n", encoding="utf-8")
    (cfg.input_dir / "target.csv").write_text("id\n2\n", encoding="utf-8")

    doc_step = PipelineStep(
        name="doc",
        main=lambda args: 0,
        input_filename="document.csv",
        output_stem="documents",
    )
    target_step = PipelineStep(
        name="target",
        main=lambda args: 0,
        input_filename="target.csv",
        output_stem="targets",
    )

    invoked_second = False

    def fail_runner(*args):
        return StepExecutionResult(exit_code=1, executed=True, status="failed", reason="boom")

    def success_runner(*args):
        nonlocal invoked_second
        invoked_second = True
        return StepExecutionResult(exit_code=0, executed=True, status="success")

    prepared = [
        PreparedPipelineStep(step=doc_step, invoke=fail_runner),
        PreparedPipelineStep(step=target_step, invoke=success_runner),
    ]

    results = list(execute_workflow(cfg, prepared))

    assert len(results) == 1
    assert results[0].result.exit_code == 1
    assert invoked_second is False


@pytest.mark.unit
def test_execute_workflow__propagates_exceptions(tmp_path: Path) -> None:
    cfg = _make_config(tmp_path)
    (cfg.input_dir / "document.csv").write_text("id\n1\n", encoding="utf-8")
    (cfg.input_dir / "target.csv").write_text("id\n2\n", encoding="utf-8")

    doc_step = PipelineStep(
        name="doc",
        main=lambda args: 0,
        input_filename="document.csv",
        output_stem="documents",
    )
    target_step = PipelineStep(
        name="target",
        main=lambda args: 0,
        input_filename="target.csv",
        output_stem="targets",
    )

    invoked_second = False

    def error_runner(*args):
        raise RuntimeError("failure")

    def second_runner(*args):
        nonlocal invoked_second
        invoked_second = True
        return StepExecutionResult(exit_code=0, executed=True, status="success")

    prepared = [
        PreparedPipelineStep(step=doc_step, invoke=error_runner),
        PreparedPipelineStep(step=target_step, invoke=second_runner),
    ]

    with pytest.raises(RuntimeError):
        list(execute_workflow(cfg, prepared))

    assert invoked_second is False
