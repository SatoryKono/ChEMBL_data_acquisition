from __future__ import annotations

import argparse
import io
import json
from dataclasses import replace
from pathlib import Path
from typing import Callable

import pandas as pd
import pytest

from scripts import get_data


def _build_stub_step(
    *,
    required_columns: list[str],
    key_column: str,
    on_execute: Callable[[pd.DataFrame, Path], int],
) -> Callable[[list[str]], int]:
    def _main(argv: list[str]) -> int:
        parser = argparse.ArgumentParser(prog="stub")
        parser.add_argument("--config", required=True)
        parser.add_argument("--input", required=True)
        parser.add_argument("--final-out", required=True)
        parser.add_argument("--log-level")
        parser.add_argument("--limit", type=int, default=None)
        parser.add_argument("--force", action="store_true")
        parser.add_argument("--skip-existing", action="store_true")
        args, _ = parser.parse_known_args(argv)
        frame = pd.read_csv(Path(args.input))
        missing = [col for col in required_columns if col not in frame.columns]
        if missing:
            get_data._LOGGER.error("schema_mismatch", missing=missing)
            return 1
        frame = frame[required_columns].copy()
        frame[key_column] = frame[key_column].astype("string").str.strip()
        duplicates = frame[key_column].duplicated()
        if duplicates.any():
            get_data._LOGGER.warning(
                "duplicates_detected",
                count=int(duplicates.sum()),
            )
            frame = frame.loc[~duplicates].copy()
        return on_execute(frame, Path(args.final_out))

    return _main


def _prepare_environment(tmp_path: Path) -> get_data.PipelineRunConfig:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
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


def _write_input(cfg: get_data.PipelineRunConfig, name: str, frame: pd.DataFrame) -> Path:
    path = cfg.input_path(name)
    frame.to_csv(path, index=False)
    return path


@pytest.mark.integration
def test_pipeline_subset__schema_and_duplicates(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = _prepare_environment(tmp_path)
    frame = pd.DataFrame(
        [
            {"document_chembl_id": "CHEMBL1", "title": "One", "pubmed_id": "1"},
            {"document_chembl_id": "CHEMBL1", "title": "Duplicate", "pubmed_id": "1"},
            {"document_chembl_id": "CHEMBL2", "title": "Two", "pubmed_id": "2"},
        ]
    )
    _write_input(cfg, "document", frame)

    output_payload: list[pd.DataFrame] = []

    def _on_execute(rows: pd.DataFrame, destination: Path) -> int:
        output_payload.append(rows)
        destination.parent.mkdir(parents=True, exist_ok=True)
        rows.to_csv(destination, index=False)
        return 0

    step = get_data.PipelineStep(
        "document",
        _build_stub_step(
            required_columns=["document_chembl_id", "title", "pubmed_id"],
            key_column="document_chembl_id",
            on_execute=_on_execute,
        ),
        None,
    )

    stream = io.StringIO()
    logger = get_data.Logger(get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration"))
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", (step,), raising=False)
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    status = get_data.run_pipeline(cfg)
    assert status == 0
    assert output_payload, "expected pipeline to write output"
    output_frame = output_payload[0]
    assert list(output_frame["document_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
    logs = [json.loads(line) for line in stream.getvalue().splitlines() if line.strip()]
    assert any(record.get("event") == "duplicates_detected" for record in logs)

    malformed = frame.drop(columns=["title"])
    _write_input(cfg, "document", malformed)
    stream.truncate(0)
    stream.seek(0)

    status_malformed = get_data.run_pipeline(cfg)
    assert status_malformed == 1
    logs = [json.loads(line) for line in stream.getvalue().splitlines() if line.strip()]
    assert any(record.get("event") == "schema_mismatch" for record in logs)


@pytest.mark.integration
def test_pipeline_subset__skip_existing_and_force(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = _prepare_environment(tmp_path)
    _write_input(
        cfg,
        "target",
        pd.DataFrame(
            [
                {"target_chembl_id": "T1", "name": "Alpha", "organism": "Human"},
                {"target_chembl_id": "T2", "name": "Beta", "organism": "Mouse"},
            ]
        ),
    )

    executions: list[int] = []

    def _on_execute(rows: pd.DataFrame, destination: Path) -> int:
        executions.append(len(rows))
        destination.write_text("target_chembl_id\nT1\nT2\n", encoding="utf-8")
        return 0

    step = get_data.PipelineStep(
        "target",
        _build_stub_step(
            required_columns=["target_chembl_id", "name", "organism"],
            key_column="target_chembl_id",
            on_execute=_on_execute,
        ),
        None,
    )

    stream = io.StringIO()
    logger = get_data.Logger(get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration"))
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", (step,), raising=False)
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    status_first = get_data.run_pipeline(cfg)
    assert status_first == 0
    assert executions == [2]

    cfg_skip = replace(cfg, skip_existing=True)
    status_skip = get_data.run_pipeline(cfg_skip)
    assert status_skip == 0
    assert executions == [2]
    logs = [json.loads(line) for line in stream.getvalue().splitlines() if line.strip()]
    assert any(record.get("event") == "step_skipped_existing" for record in logs)

    cfg_force = replace(cfg, skip_existing=True, force=True)
    status_force = get_data.run_pipeline(cfg_force)
    assert status_force == 0
    assert executions == [2, 2]


@pytest.mark.integration
def test_pipeline_subset__retry_after_failure(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = _prepare_environment(tmp_path)
    _write_input(
        cfg,
        "assay",
        pd.DataFrame(
            [
                {
                    "assay_chembl_id": "A1",
                    "target_chembl_id": "T1",
                    "document_chembl_id": "D1",
                    "description": "First",
                }
            ]
        ),
    )

    attempts = {"count": 0}

    def _on_execute(rows: pd.DataFrame, destination: Path) -> int:
        attempts["count"] += 1
        if attempts["count"] == 1:
            destination.parent.mkdir(parents=True, exist_ok=True)
            tmp_path = destination.with_suffix(".tmp")
            tmp_path.write_text("partial\n", encoding="utf-8")
            return 1
        destination.write_text("assay_chembl_id\nA1\n", encoding="utf-8")
        return 0

    step = get_data.PipelineStep(
        "assay",
        _build_stub_step(
            required_columns=[
                "assay_chembl_id",
                "target_chembl_id",
                "document_chembl_id",
                "description",
            ],
            key_column="assay_chembl_id",
            on_execute=_on_execute,
        ),
        None,
    )

    stream = io.StringIO()
    logger = get_data.Logger(get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration"))
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", (step,), raising=False)
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    first_status = get_data.run_pipeline(cfg)
    assert first_status == 1
    working = get_data._temporary_output_path(step.expected_output(cfg))
    assert not working.exists()
    sentinel = get_data._failure_sentinel_path(step.expected_output(cfg))
    assert sentinel.exists()

    sentinel.unlink()
    second_status = get_data.run_pipeline(cfg)
    assert second_status == 0
    final_output = step.expected_output(cfg)
    assert final_output.exists()
    assert attempts["count"] == 2
