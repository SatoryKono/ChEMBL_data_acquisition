"""Unit tests for orchestration artefact helpers."""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from library.cli.commands import get_data
from library.common.logging_setup import LoggerConfig, configure_logger
from library.orchestration import artifacts


@pytest.mark.unit
@pytest.mark.parametrize(
    "filename",
    [
        "output.targets_20240101_failure_cases.csv",
        "output.targets_20240101_chembl.csv",
        "output.targets_20240101_uniprot.csv",
        "output.documents_20240101_pubmed.csv",
        "output.targets_20240101_normalized.csv",
        "output.targets_raw.csv",
        "output.targets.raw.parquet",
        "output.targets_20240101.quality.json",
        "output.targets_20240101.postprocess.report.json",
        "output_postprocessed.targets.csv",
    ],
)
def test_is_diagnostic_sidecar__legacy_suffixes(filename: str) -> None:
    assert artifacts.is_diagnostic_sidecar(Path(filename)) is True


@pytest.mark.unit
@pytest.mark.parametrize(
    "filename",
    [
        "output.targets_20240101.csv_chembl_20240101.csv",
        "output.targets_20240101.csv_chembl_20240101_data_correlation_report_table.csv",
        "output.targets_20240101.csv_chembl_20240101_quality_report_table.csv",
    ],
)
def test_is_diagnostic_sidecar__intermediate_standard_outputs(filename: str) -> None:
    assert artifacts.is_diagnostic_sidecar(Path(filename)) is True


@pytest.mark.unit
@pytest.mark.parametrize(
    "filename",
    [
        "output.targets_20240101.csv",
        "output.targets_20240101_quality_report_table.csv",
        "output.targets_20240101_data_correlation_report_table.csv",
    ],
)
def test_is_diagnostic_sidecar__canonical_outputs(filename: str) -> None:
    assert artifacts.is_diagnostic_sidecar(Path(filename)) is False


@pytest.mark.unit
def test_run_postprocess_hook__missing_input(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    step = next(step for step in get_data.DEFAULT_PIPELINE_STEPS if step.name == "document")
    final_output = tmp_path / "output.documents_20240101.csv"
    assert not final_output.exists()

    calls: list[tuple[str, dict[str, object]]] = []

    def _record(event: str, **kwargs: object) -> None:
        calls.append((event, kwargs))

    logger = configure_logger(LoggerConfig(level="INFO", stream=io.StringIO()))
    monkeypatch.setattr(logger, "warning", _record)

    result = artifacts.run_postprocess_hook(
        step,
        final_output,
        table="documents",
        logger=logger,
    )

    assert result is None
    assert calls == [
        (
            "postprocess_input_missing",
            {
                "step": step.name,
                "table": "documents",
                "input": str(final_output),
            },
        )
    ]


@pytest.mark.unit
def test_finalize_outputs_removes_diagnostics(tmp_path: Path) -> None:
    working_dir = tmp_path / "work"
    final_dir = tmp_path / "final"
    working_output = working_dir / "output.targets_20240101.csv"
    final_output = final_dir / working_output.name
    sentinel = final_output.with_name(f"{final_output.name}.failed")

    working_output.parent.mkdir(parents=True, exist_ok=True)
    final_dir.mkdir(parents=True, exist_ok=True)
    working_output.write_text("content")
    diagnostic_sidecar = working_output.with_name("output.targets_20240101_failure_cases.csv")
    diagnostic_sidecar.write_text("diag")
    sentinel.touch()

    logger = configure_logger(LoggerConfig(level="INFO", stream=io.StringIO()))

    artifacts.finalize_outputs(
        final_output,
        working_output,
        sentinel,
        diagnostics_enabled=False,
        logger=logger,
        include_patterns=("*.csv",),
    )

    assert final_output.exists()
    assert final_output.read_text() == "content"
    assert not diagnostic_sidecar.exists()
    assert not sentinel.exists()
