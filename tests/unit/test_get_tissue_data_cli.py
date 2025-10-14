"""Unit tests for the :mod:`scripts.get_tissue_data` CLI helpers."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pytest
from scripts import get_tissue_data

from library.config import Config
from library.pipelines.tissue.pipeline import TissuePipelineResult


class _StubLogger:
    """Collects structured logging events for assertions."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def bind(self, **_: object) -> _StubLogger:  # pragma: no cover - interface parity
        return self

    def info(self, event: str, **data: object) -> None:
        self.events.append(("info", event, dict(data)))

    def warning(
        self, event: str, **data: object
    ) -> None:  # pragma: no cover - defensive
        self.events.append(("warning", event, dict(data)))

    def error(self, event: str, **data: object) -> None:  # pragma: no cover - defensive
        self.events.append(("error", event, dict(data)))


@pytest.mark.unit
def test_get_tissue_data_main__limit_zero_skips_pipeline(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """``main()`` should short-circuit when ``--limit 0`` is provided."""

    stub_logger = _StubLogger()
    monkeypatch.setattr(get_tissue_data, "logger", stub_logger)

    def _unexpected_run_cli_command(**_: object) -> None:
        raise AssertionError("run_cli_command should not be invoked when limit is zero")

    monkeypatch.setattr(get_tissue_data, "run_cli_command", _unexpected_run_cli_command)
    monkeypatch.setattr(
        get_tissue_data, "configure_logger", lambda *_args, **_kwargs: None
    )

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("tissue_chembl_id\nCHEMBLT0\n", encoding="utf-8")
    final_out = tmp_path / "output.csv"

    exit_code = get_tissue_data.main(
        [
            "--input",
            str(input_csv),
            "--final-out",
            str(final_out),
            "--limit",
            "0",
        ]
    )

    assert exit_code == 0
    assert ("info", "pipeline_skip_limit", {"limit": 0}) in stub_logger.events
    assert final_out.exists() is False


@pytest.mark.unit
def test_run__persists_metadata_sidecar(
    cfg: Config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("tissue_chembl_id\nCHEMBLT0\n", encoding="utf-8")
    output_csv = tmp_path / "output.tissue_20200101.csv"
    output_csv.write_text("tissue_chembl_id\nCHEMBLT0\n", encoding="utf-8")

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        skip_existing=False,
        force=False,
        offset=0,
        limit=None,
        mode="chembl",
        batch_size=5,
        timeout=None,
        column="tissue_chembl_id",
    )

    result = TissuePipelineResult(
        exit_code=0,
        records=1,
        duration=0.2,
        output_path=output_csv,
        failure_path=None,
        failures=0,
        missing_ids=(),
        written=True,
    )

    monkeypatch.setattr(get_tissue_data, "run_tissue_pipeline", lambda *_: result)
    monkeypatch.setattr(
        get_tissue_data,
        "generate_correlation_report",
        lambda *_, **__: pd.DataFrame(),
    )
    monkeypatch.setattr(
        get_tissue_data,
        "generate_qc_report",
        lambda *_, **__: pd.DataFrame(),
    )

    metadata_calls: list[tuple[tuple[object, ...], dict[str, object]]] = []
    original_save_metadata = get_tissue_data.io.save_metadata

    def _tracking_save_metadata(*args: object, **kwargs: object):
        metadata_calls.append((args, kwargs))
        return original_save_metadata(*args, **kwargs)

    monkeypatch.setattr(get_tissue_data.io, "save_metadata", _tracking_save_metadata)

    exit_code = get_tissue_data.run(cfg, args)

    assert exit_code == 0
    assert metadata_calls, "save_metadata should be invoked"

    _, call_kwargs = metadata_calls[-1]
    assert call_kwargs["qc_summary"] == {"rows": 1}

    artifacts = call_kwargs["artifacts"]
    assert isinstance(artifacts, list)
    dataset_artifact = args.output_csv
    expected_artifact_names = {
        dataset_artifact.name,
        f"{dataset_artifact.stem}_quality_report_table.csv",
        f"{dataset_artifact.stem}_data_correlation_report_table.csv",
    }
    assert {Path(path).name for path in artifacts} == expected_artifact_names

    table_name = str(call_kwargs["table_name"])
    date_tag = str(call_kwargs["date_tag"])
    meta_path = dataset_artifact.parent / f"output.{table_name}_{date_tag}.meta.yaml"
    assert meta_path.exists()
