from __future__ import annotations

import argparse
import io
from dataclasses import replace
from pathlib import Path
from typing import Callable

import pandas as pd
import pytest

from library.config import Config
from scripts import get_data, get_target_data
from tests.helpers.logs import parse_log_lines


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
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration")
    )
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", (step,), raising=False)
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    status = get_data.run_pipeline(cfg)
    assert status == 0
    assert output_payload, "expected pipeline to write output"
    output_frame = output_payload[0]
    assert list(output_frame["document_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
    logs = parse_log_lines(stream.getvalue())
    assert any(record.get("event") == "duplicates_detected" for record in logs)

    malformed = frame.drop(columns=["title"])
    _write_input(cfg, "document", malformed)
    stream.truncate(0)
    stream.seek(0)

    status_malformed = get_data.run_pipeline(cfg)
    assert status_malformed == 1
    logs = parse_log_lines(stream.getvalue())
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
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration")
    )
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", (step,), raising=False)
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    status_first = get_data.run_pipeline(cfg)
    assert status_first == 0
    assert executions == [2]

    cfg_skip = replace(cfg, skip_existing=True)
    status_skip = get_data.run_pipeline(cfg_skip)
    assert status_skip == 0
    assert executions == [2]
    logs = parse_log_lines(stream.getvalue())
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
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration")
    )
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


@pytest.mark.integration
def test_pipeline_subset__target_postprocess_sidecars(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = _prepare_environment(tmp_path)
    _write_input(
        cfg,
        "target",
        pd.DataFrame(
            [
                {
                    "target_chembl_id": "CHEMBL1",
                    "target_name": "Alpha",
                    "organism": "Homo sapiens",
                }
            ],
            dtype="string",
        ),
    )

    call_order: list[str] = []

    def _sidecar_path(source: Path, prefix: str) -> Path:
        canonical = get_target_data._normalise_target_export_name(source).lstrip(".")
        destination = source.with_name(f"{prefix}{canonical}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
        return destination

    def _fake_organism(source: Path, *, cfg: Config) -> Path:
        call_order.append("organism")
        return _sidecar_path(source, "organism.")

    def _fake_isoform(
        source: Path,
        *,
        cfg: Config,
        context: get_target_data.IsoformPostprocessContext | None = None,
        ambiguous_classifications: int | None = None,
    ) -> Path:
        call_order.append("isoform")
        return _sidecar_path(source, "isoform.")

    def _fake_names(source: Path, *, cfg: Config) -> Path:
        call_order.append("names")
        return _sidecar_path(source, "name.")

    def _fake_iuphar(source: Path, *, verbose: bool = True) -> Path:
        call_order.append("iuphar")
        return _sidecar_path(source, "IUPHAR.")

    def _stub_run_all(cfg_obj: Config, args: argparse.Namespace) -> int:
        working_output = Path(args.final_out)
        working_output.parent.mkdir(parents=True, exist_ok=True)
        working_output.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
        get_target_data._postprocess_target_exports(working_output, cfg=cfg_obj)
        return 0

    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration")
    )

    step = get_data.PipelineStep("target", get_target_data.main, "all")

    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", (step,), raising=False)
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)
    monkeypatch.setattr(get_target_data, "run_all", _stub_run_all)
    monkeypatch.setattr(
        get_target_data,
        "_postprocess_organism_export",
        _fake_organism,
    )
    monkeypatch.setattr(get_target_data, "_postprocess_isoform_export", _fake_isoform)
    monkeypatch.setattr(get_target_data, "_postprocess_names_export", _fake_names)
    monkeypatch.setattr(get_target_data, "_postprocess_iuphar_export", _fake_iuphar)

    status = get_data.run_pipeline(cfg)

    assert status == 0
    assert call_order == ["organism", "isoform", "names", "iuphar"]

    final_output = step.expected_output(cfg)
    assert final_output.exists()
    sidecars = {
        "organism": final_output.with_name(f"organism.{final_output.name}"),
        "isoform": final_output.with_name(f"isoform.{final_output.name}"),
        "names": final_output.with_name(f"name.{final_output.name}"),
        "iuphar": final_output.with_name(f"IUPHAR.{final_output.name}"),
    }
    for path in sidecars.values():
        assert path.exists()

