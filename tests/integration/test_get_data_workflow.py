from __future__ import annotations

import argparse
import importlib
import io
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from pydantic import BaseModel, ValidationError

from library.config import Config, ConfigLoaderError
from library.pipelines.common import PipelineRunResult
from scripts import get_target_data
from tests.helpers import ASSAY_ENRICHMENT_MIN_RATIO
from tests.helpers.logs import parse_log_lines
from tests.helpers.manifests import list_manifest_files, load_latest_manifest


get_data = importlib.import_module("library.cli.commands.get_data")


def _build_stub_api(
    *,
    required_columns: list[str],
    key_column: str,
    on_execute: Callable[[pd.DataFrame, Path], int],
) -> get_data.PipelineApi:
    def _builder(
        cfg: get_data.PipelineRunConfig, input_path: Path, output_path: Path
    ) -> SimpleNamespace:
        return SimpleNamespace(input_csv=input_path, output_csv=output_path)

    def _runner(config: Config, options: SimpleNamespace) -> PipelineRunResult:
        frame = pd.read_csv(Path(options.input_csv))
        missing = [col for col in required_columns if col not in frame.columns]
        if missing:
            get_data._LOGGER.error("schema_mismatch", missing=missing)
            return PipelineRunResult(
                exit_code=1,
                output_path=Path(options.output_csv),
                executed=True,
                reason="schema_mismatch",
                written=False,
            )
        frame = frame[required_columns].copy()
        frame[key_column] = frame[key_column].astype("string").str.strip()
        duplicates = frame[key_column].duplicated()
        if duplicates.any():
            get_data._LOGGER.warning(
                "duplicates_detected",
                count=int(duplicates.sum()),
            )
            frame = frame.loc[~duplicates].copy()
        destination = Path(options.output_csv)
        destination.parent.mkdir(parents=True, exist_ok=True)
        result_code = on_execute(frame, destination)
        return PipelineRunResult(
            exit_code=result_code,
            output_path=destination,
            executed=True,
            reason=None if result_code == 0 else "pipeline_failed",
            written=result_code == 0,
        )

    return get_data.PipelineApi(_builder, _runner)


def _prepare_environment(tmp_path: Path) -> get_data.PipelineRunConfig:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    config_path = base_path / "config.yaml"
    config_path.write_text("io:\n  csv_sep: ','\n", encoding="utf-8")
    input_files = dict(get_data.DEFAULT_INPUT_FILES)
    output_stems = dict(get_data.DEFAULT_OUTPUT_STEMS)
    subcommands = {
        step.name: step.subcommand for step in get_data.DEFAULT_PIPELINE_STEPS
    }
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
        input_files=input_files,
        output_stems=output_stems,
        subcommands=subcommands,
    )


def _write_input(
    cfg: get_data.PipelineRunConfig, name: str, frame: pd.DataFrame
) -> Path:
    path = cfg.input_path(name)
    frame.to_csv(path, index=False)
    return path


def _load_manifest(cfg: get_data.PipelineRunConfig) -> dict[str, object]:
    _, manifest = load_latest_manifest(cfg.base_path)
    return manifest


def _make_validation_error() -> ValidationError:
    class _DummyModel(BaseModel):
        value: int

    try:
        _DummyModel.model_validate({"value": "boom"})
    except ValidationError as exc:  # pragma: no cover - control flow
        return exc
    raise AssertionError("expected ValidationError")


@pytest.mark.integration
def test_pipeline_subset__schema_and_duplicates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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

    api = _build_stub_api(
        required_columns=["document_chembl_id", "title", "pubmed_id"],
        key_column="document_chembl_id",
        on_execute=_on_execute,
    )
    step = get_data.PipelineStep(
        name="document",
        main=lambda _: 0,
        input_filename="document.csv",
        output_stem="documents",
    )

    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration")
    )
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)
    monkeypatch.setattr(
        get_data,
        "load_config",
        lambda *args, **kwargs: Config(),
        raising=False,
    )
    monkeypatch.setattr(get_data, "_PIPELINE_APIS", {"document": api}, raising=False)

    steps = (step,)
    status = get_data.run_pipeline(cfg, steps=steps)
    assert status == 0
    assert output_payload, "expected pipeline to write output"
    output_frame = output_payload[0]
    assert list(output_frame["document_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
    logs = parse_log_lines(stream.getvalue())
    assert any(record.get("event") == "duplicates_detected" for record in logs)
    manifests_after_success = list_manifest_files(cfg.base_path)
    assert len(manifests_after_success) == 1
    manifest_success = _load_manifest(cfg)
    assert manifest_success["run"]["exit_code"] == 0
    assert manifest_success["steps"][0]["status"] == "success"
    assert manifest_success["steps"][0]["output"]["exists"] is True
    assert manifest_success["steps"][0]["output"]["checksum_sha256"]

    malformed = frame.drop(columns=["title"])
    _write_input(cfg, "document", malformed)
    stream.truncate(0)
    stream.seek(0)

    status_malformed = get_data.run_pipeline(cfg, steps=steps)
    assert status_malformed == 1
    logs = parse_log_lines(stream.getvalue())
    assert any(record.get("event") == "schema_mismatch" for record in logs)
    manifests_after_failure = list_manifest_files(cfg.base_path)
    assert len(manifests_after_failure) == 2
    manifest_failure = _load_manifest(cfg)
    assert manifest_failure["run"]["exit_code"] == 1
    assert manifest_failure["steps"][0]["status"] == "failed"
    assert manifest_failure["steps"][0]["reason"] == "schema_mismatch"
    assert manifest_failure["steps"][0]["output"]["exists"] is False


@pytest.mark.integration
def test_run_pipeline__config_loader_error_handled(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = _prepare_environment(tmp_path)
    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="config-error")
    )
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    error = ConfigLoaderError("failed to load configuration")

    def _raise_config_error(*_args: object, **_kwargs: object) -> Config:
        raise error

    monkeypatch.setattr(get_data, "load_config", _raise_config_error, raising=False)

    status = get_data.run_pipeline(cfg, steps=())

    assert status == 1
    log_text = stream.getvalue()
    assert "config_load_failed" in log_text
    assert "pipeline_done exit_code=1" in log_text
    assert list_manifest_files(cfg.base_path) == []


@pytest.mark.integration
def test_run_pipeline__validation_error_handled(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = _prepare_environment(tmp_path)
    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="validation-error")
    )
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    validation_error = _make_validation_error()

    def _raise_validation_error(*_args: object, **_kwargs: object) -> Config:
        raise validation_error

    monkeypatch.setattr(get_data, "load_config", _raise_validation_error, raising=False)

    status = get_data.run_pipeline(cfg, steps=())

    assert status == 1
    log_text = stream.getvalue()
    assert "config_load_failed" in log_text
    assert "pipeline_done exit_code=1" in log_text
    assert list_manifest_files(cfg.base_path) == []


@pytest.mark.integration
def test_pipeline_subset__skip_existing_and_force(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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

    api = _build_stub_api(
        required_columns=["target_chembl_id", "name", "organism"],
        key_column="target_chembl_id",
        on_execute=_on_execute,
    )
    step = get_data.PipelineStep(
        name="target",
        main=lambda _: 0,
        input_filename="target.csv",
        output_stem="targets",
    )

    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration")
    )
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    steps = (step,)
    monkeypatch.setattr(get_data, "_PIPELINE_APIS", {"target": api}, raising=False)
    status_first = get_data.run_pipeline(cfg, steps=steps)
    assert status_first == 0
    assert executions == [2]
    manifest_first = _load_manifest(cfg)
    assert manifest_first["steps"][0]["status"] == "success"
    assert manifest_first["steps"][0]["output"]["exists"] is True

    cfg_skip = replace(cfg, skip_existing=True)
    status_skip = get_data.run_pipeline(cfg_skip, steps=steps)
    assert status_skip == 0
    assert executions == [2]
    logs = parse_log_lines(stream.getvalue())
    assert any(record.get("event") == "step_skipped_existing" for record in logs)
    manifest_skip = _load_manifest(cfg_skip)
    skip_entry = manifest_skip["steps"][0]
    assert skip_entry["status"] == "skipped"
    assert skip_entry["reason"] == "skip_existing"
    assert skip_entry["executed"] is False
    assert skip_entry["output"]["exists"] is True

    cfg_force = replace(cfg, skip_existing=True, force=True)
    status_force = get_data.run_pipeline(cfg_force, steps=steps)
    assert status_force == 0
    assert executions == [2, 2]
    manifest_force = _load_manifest(cfg_force)
    assert manifest_force["steps"][0]["status"] == "success"
    assert manifest_force["steps"][0]["output"]["exists"] is True


@pytest.mark.integration
def test_run_pipeline__target_override_invokes_target_branch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    config_path = base_path / "config.yaml"
    config_path.write_text("io:\n  csv_sep: ','\n", encoding="utf-8")

    args = SimpleNamespace(
        base_path=base_path,
        input_dir=Path("input"),
        output_dir=Path("output"),
        config=config_path,
        date_prefix="20200101",
        log_level="INFO",
        verbose=False,
        limit=None,
        force=False,
        skip_existing=False,
        dry_run=False,
        pipeline_registry=None,
        override_input=[],
        override_output_stem=[],
        override_subcommand=["target=uniprot"],
    )

    steps = get_data._resolve_pipeline_steps(args)
    cfg = get_data._prepare_config(args, steps)
    assert cfg.subcommand_for("target") == "uniprot"

    target_input = cfg.input_path("target")
    pd.DataFrame(
        [
            {"target_chembl_id": "CHEMBLT1", "name": "Target", "organism": "Human"},
        ]
    ).to_csv(target_input, index=False)

    target_step = next(step for step in steps if step.name == "target")

    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration")
    )
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)
    monkeypatch.setattr(
        get_data,
        "load_config",
        lambda *args, **kwargs: Config(),
        raising=False,
    )

    commands: list[str] = []

    def _fake_run_target_pipeline(config: Config, options: object) -> PipelineRunResult:
        command = options.command
        commands.append(command)
        assert command == "uniprot"
        destination = Path(options.output_csv)
        destination.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            [
                {
                    "target_chembl_id": "CHEMBLT1",
                    "name": "Target",
                    "organism": "Human",
                }
            ]
        ).to_csv(destination, index=False)
        return PipelineRunResult(
            exit_code=0,
            output_path=destination,
            executed=True,
            reason=None,
            written=True,
        )

    monkeypatch.setattr(
        get_data,
        "_PIPELINE_APIS",
        {
            "target": get_data.PipelineApi(
                get_data._build_target_options, _fake_run_target_pipeline
            )
        },
        raising=False,
    )

    status = get_data.run_pipeline(cfg, steps=(target_step,))

    assert status == 0
    assert commands == ["uniprot"]
    final_output = cfg.output_path("target")
    assert final_output.exists()


@pytest.mark.integration
def test_run_pipeline__document_override_invokes_selected_mode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    config_path = base_path / "config.yaml"
    config_path.write_text("io:\n  csv_sep: ','\n", encoding="utf-8")

    args = SimpleNamespace(
        base_path=base_path,
        input_dir=Path("input"),
        output_dir=Path("output"),
        config=config_path,
        date_prefix="20200101",
        log_level="INFO",
        verbose=False,
        limit=None,
        force=False,
        skip_existing=False,
        dry_run=False,
        pipeline_registry=None,
        override_input=[],
        override_output_stem=[],
        override_subcommand=["document=chembl"],
    )

    steps = get_data._resolve_pipeline_steps(args)
    document_step = next(step for step in steps if step.name == "document")
    cfg = get_data._prepare_config(args, (document_step,))
    assert cfg.subcommand_for("document") == "chembl"

    document_input = cfg.input_path("document")
    pd.DataFrame(
        [
            {
                "document_chembl_id": "CHEMBLDOC1",
                "title": "Document",
            }
        ]
    ).to_csv(document_input, index=False)

    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration_doc")
    )
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)
    monkeypatch.setattr(
        get_data,
        "load_config",
        lambda *args, **kwargs: Config(),
        raising=False,
    )

    modes: list[str] = []

    def _fake_run_document_pipeline(config: Config, options: object) -> PipelineRunResult:
        mode = getattr(options, "mode")
        modes.append(mode)
        assert mode == "chembl"
        destination = Path(options.output_csv)
        destination.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            [
                {
                    "document_chembl_id": "CHEMBLDOC1",
                    "title": "Document",
                }
            ]
        ).to_csv(destination, index=False)
        return PipelineRunResult(
            exit_code=0,
            output_path=destination,
            executed=True,
            reason=None,
            written=True,
        )

    monkeypatch.setattr(
        get_data,
        "_PIPELINE_APIS",
        {
            "document": get_data.PipelineApi(
                get_data._build_document_options, _fake_run_document_pipeline
            )
        },
        raising=False,
    )

    status = get_data.run_pipeline(cfg, steps=(document_step,))

    assert status == 0
    assert modes == ["chembl"]
    final_output = cfg.output_path("document")
    assert final_output.exists()


@pytest.mark.integration
def test_pipeline_subset__testitem_skip_existing_avoids_parent_warm(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = _prepare_environment(tmp_path)
    _write_input(
        cfg,
        "testitem",
        pd.DataFrame(
            [
                {"testitem_chembl_id": "CHEMBLT1", "compound_key": "A"},
                {"testitem_chembl_id": "CHEMBLT2", "compound_key": "B"},
            ],
            dtype="string",
        ),
    )

    executions: list[int] = []

    def _on_execute(rows: pd.DataFrame, destination: Path) -> int:
        executions.append(len(rows))
        destination.parent.mkdir(parents=True, exist_ok=True)
        rows.to_csv(destination, index=False)
        return 0

    api = _build_stub_api(
        required_columns=["testitem_chembl_id", "compound_key"],
        key_column="testitem_chembl_id",
        on_execute=_on_execute,
    )

    step = get_data.PipelineStep(
        name="testitem",
        main=lambda _: 0,
        input_filename="testitem.csv",
        output_stem="testitems",
    )

    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration")
    )

    warm_calls: list[tuple[get_data.PipelineRunConfig, object]] = []

    def _record_warm(
        current_cfg: get_data.PipelineRunConfig, base_config: object
    ) -> None:
        warm_calls.append((current_cfg, base_config))

    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)
    monkeypatch.setattr(get_data, "_PIPELINE_APIS", {"testitem": api}, raising=False)
    monkeypatch.setattr(get_data, "_warm_parent_catalog", _record_warm, raising=False)

    status_first = get_data.run_pipeline(cfg, steps=(step,))
    assert status_first == 0
    assert executions == [2]
    assert len(warm_calls) == 1

    final_output = step.expected_output(cfg)
    assert final_output.exists()

    cfg_skip = replace(cfg, skip_existing=True)
    status_second = get_data.run_pipeline(cfg_skip, steps=(step,))
    assert status_second == 0
    assert executions == [2]
    assert len(warm_calls) == 1

    logs = parse_log_lines(stream.getvalue())
    assert any(
        record.get("event") == "step_skipped_existing"
        and record.get("data", {}).get("step") == "testitem"
        for record in logs
    )

    manifest_skip = _load_manifest(cfg_skip)
    skip_entry = manifest_skip["steps"][0]
    assert skip_entry["status"] == "skipped"
    assert skip_entry["reason"] == "skip_existing"
    assert skip_entry["executed"] is False


@pytest.mark.integration
def test_pipeline_subset__retry_after_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = _prepare_environment(tmp_path)
    _write_input(
        cfg,
        "assay",
        pd.DataFrame(
            [
                {
                    "assay_chembl_id": "CHEMBLA1",
                    "target_chembl_id": "CHEMBLT1",
                    "document_chembl_id": "CHEMBL123",
                    "description": "Binding assay",
                }
            ]
        ),
    )

    attempts = {"count": 0}
    dictionary_path = (
        Path(__file__).resolve().parents[1] / "resources" / "pipeline_inputs" / "assay_dictionary.csv"
    )

    def _on_execute(rows: pd.DataFrame, destination: Path) -> int:
        attempts["count"] += 1
        if attempts["count"] == 1:
            destination.parent.mkdir(parents=True, exist_ok=True)
            tmp_path = destination.with_suffix(".tmp")
            tmp_path.write_text("partial\n", encoding="utf-8")
            return 1
        dictionary = pd.read_csv(dictionary_path)
        dictionary["assay_chembl_id"] = dictionary["assay_chembl_id"].astype("string")
        enriched = rows.merge(dictionary, on="assay_chembl_id", how="left")
        enriched["description"] = enriched["description"].astype("string").str.strip()
        enriched["description_length"] = (
            enriched["description"].str.len().astype("Int64")
        )
        enriched["year"] = pd.to_numeric(enriched["year"], errors="coerce").astype(
            "Int64"
        )
        quality_columns = ["assay_strain", "assay_group", "year", "accession"]
        completeness = 1.0 - enriched[quality_columns].isna().mean()
        assert (
            completeness >= ASSAY_ENRICHMENT_MIN_RATIO
        ).all(), completeness.to_dict()
        columns = [
            "assay_chembl_id",
            "target_chembl_id",
            "document_chembl_id",
            "description",
            "description_length",
            "assay_strain",
            "assay_group",
            "year",
            "accession",
        ]
        enriched.to_csv(destination, index=False, columns=columns)
        return 0

    api = _build_stub_api(
        required_columns=[
            "assay_chembl_id",
            "target_chembl_id",
            "document_chembl_id",
            "description",
        ],
        key_column="assay_chembl_id",
        on_execute=_on_execute,
    )
    step = get_data.PipelineStep(
        name="assay",
        main=lambda _: 0,
        input_filename="assay.csv",
        output_stem="assays",
    )

    stream = io.StringIO()
    logger = get_data.configure_logger(
        get_data.LoggerConfig(level="DEBUG", stream=stream, run_id="integration")
    )
    monkeypatch.setattr(get_data, "_LOGGER", logger, raising=False)

    monkeypatch.setattr(get_data, "_PIPELINE_APIS", {"assay": api}, raising=False)
    steps = (step,)
    first_status = get_data.run_pipeline(cfg, steps=steps)
    assert first_status == 1
    working = get_data._temporary_output_path(step.expected_output(cfg))
    assert not working.exists()
    sentinel = get_data._failure_sentinel_path(step.expected_output(cfg))
    assert sentinel.exists()
    manifest_failure = _load_manifest(cfg)
    failure_entry = manifest_failure["steps"][0]
    assert failure_entry["status"] == "failed"
    assert failure_entry["reason"] == "pipeline_failed"
    assert failure_entry["output"]["exists"] is False

    sentinel.unlink()
    second_status = get_data.run_pipeline(cfg, steps=steps)
    assert second_status == 0
    final_output = step.expected_output(cfg)
    assert final_output.exists()
    assert attempts["count"] == 2
    manifest_success = _load_manifest(cfg)
    assert manifest_success["steps"][0]["status"] == "success"
    assert manifest_success["steps"][0]["output"]["exists"] is True
    result = pd.read_csv(final_output)
    expected_columns = [
        "assay_chembl_id",
        "target_chembl_id",
        "document_chembl_id",
        "description",
        "description_length",
        "assay_strain",
        "assay_group",
        "year",
        "accession",
    ]
    assert list(result.columns) == expected_columns
    quality_columns = ["assay_strain", "assay_group", "year", "accession"]
    completeness = 1.0 - result[quality_columns].isna().mean()
    assert (completeness >= ASSAY_ENRICHMENT_MIN_RATIO).all(), completeness.to_dict()


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

    step = get_data.PipelineStep(
        name="target",
        main=get_target_data.main,
        input_filename="target.csv",
        output_stem="targets",
        subcommand="all",
    )

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

    status = get_data.run_pipeline(cfg, steps=(step,))

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
    manifest = _load_manifest(cfg)
    step_entry = manifest["steps"][0]
    assert step_entry["status"] == "success"
    assert step_entry["output"]["exists"] is True
    recorded_sidecars = {item["path"]: item for item in step_entry["sidecars"]}
    assert len(recorded_sidecars) == 4
    for path in sidecars.values():
        meta = recorded_sidecars[str(path)]
        assert meta["exists"] is True
        assert meta["checksum_sha256"]
