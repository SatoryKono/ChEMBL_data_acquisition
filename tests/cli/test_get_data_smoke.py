from __future__ import annotations

import builtins
import importlib
import io
from collections.abc import Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library.io.output_writer import save_standard_outputs
from library.pipelines.registry import PipelineStep


_HTTP_PAYLOADS: dict[str, list[dict[str, Any]]] = {
    "documents": [
        {"identifier": "DOC-001", "title": "Alpha", "topic": "biology"},
        {"identifier": "DOC-002", "title": "Beta", "topic": None},
    ],
    "targets": [
        {"identifier": "TGT-001", "symbol": "T1", "organism": "human"},
        {"identifier": "TGT-002", "symbol": "T2", "organism": "mouse"},
    ],
    "assays": [
        {"assay_id": "A-1", "name": "Binding", "readout": "ic50"},
        {"assay_id": "A-2", "name": "Reporter", "readout": "ec50"},
    ],
    "testitems": [
        {"sample_id": "S-1", "batch": "B1", "purity": "0.9"},
        {"sample_id": "S-2", "batch": "B2", "purity": None},
    ],
    "activities": [
        {"activity_id": "ACT-1", "target": "TGT-001", "value": "42"},
        {"activity_id": "ACT-2", "target": "TGT-002", "value": "37"},
    ],
}


@dataclass(frozen=True)
class _StepSpec:
    name: str
    endpoint: str
    columns: tuple[str, ...]


def _make_step(spec: _StepSpec) -> PipelineStep:
    def _main(argv: Sequence[str] | None = None) -> int:
        return 0

    return PipelineStep(
        name=spec.name,
        main=_main,
        input_filename=f"{spec.endpoint}.csv",
        output_stem=spec.name,
    )


def _patch_pipeline_defaults(
    get_data, monkeypatch: pytest.MonkeyPatch, steps: tuple[PipelineStep, ...]
) -> None:
    inputs = {step.name: step.input_filename for step in steps}
    outputs = {step.name: step.output_stem for step in steps}
    subcommands = {step.name: step.subcommand for step in steps}

    monkeypatch.setattr(get_data, "DEFAULT_PIPELINE_STEPS", steps)
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", steps)
    monkeypatch.setattr(
        get_data,
        "DEFAULT_INPUT_FILES",
        get_data.PipelineInputFiles.from_mapping(inputs),
    )
    monkeypatch.setattr(get_data, "_DEFAULT_INPUT_FILES", get_data.DEFAULT_INPUT_FILES)
    monkeypatch.setattr(
        get_data,
        "DEFAULT_OUTPUT_STEMS",
        get_data.PipelineOutputStems.from_mapping(outputs),
    )
    monkeypatch.setattr(get_data, "_DEFAULT_OUTPUT_STEMS", get_data.DEFAULT_OUTPUT_STEMS)
    monkeypatch.setattr(
        get_data,
        "DEFAULT_SUBCOMMANDS",
        get_data.PipelineSubcommands.from_mapping(subcommands),
    )
    monkeypatch.setattr(get_data, "_DEFAULT_SUBCOMMANDS", get_data.DEFAULT_SUBCOMMANDS)
    monkeypatch.setattr(get_data, "_PIPELINE_APIS", {})


@pytest.mark.integration
def test_get_data_smoke_creates_expected_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pandera_module = importlib.import_module("pandera")
    monkeypatch.setattr(builtins, "Column", pandera_module.Column, raising=False)
    monkeypatch.setattr(
        builtins,
        "DataFrameSchema",
        pandera_module.DataFrameSchema,
        raising=False,
    )
    get_data = importlib.import_module("library.cli.commands.get_data")
    logging_module = importlib.import_module("library.cli.logging")
    CLILoggingContext = logging_module.CLILoggingContext

    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()

    config_path = base_path / "config.yaml"
    config_path.write_text("io: {}\n", encoding="utf-8")

    specs = (
        _StepSpec("smoke_documents", "documents", ("identifier", "title", "topic")),
        _StepSpec("smoke_targets", "targets", ("identifier", "symbol", "organism")),
        _StepSpec("smoke_assays", "assays", ("assay_id", "name", "readout")),
        _StepSpec("smoke_testitems", "testitems", ("sample_id", "batch", "purity")),
        _StepSpec("smoke_activities", "activities", ("activity_id", "target", "value")),
    )
    steps = tuple(_make_step(spec) for spec in specs)

    for spec in specs:
        (input_dir / f"{spec.endpoint}.csv").write_text("identifier\nplaceholder\n", encoding="utf-8")

    _patch_pipeline_defaults(get_data, monkeypatch, steps)

    @contextmanager
    def _fake_logging(*args, **kwargs):  # noqa: ANN002, ANN003
        log_path = base_path / "logs" / "cli.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        cfg = get_data.LoggerConfig(level="INFO", run_id="smoke-run")
        yield CLILoggingContext(
            log_path=log_path,
            log_cfg=cfg,
            console_stream=io.StringIO(),
        )

    monkeypatch.setattr(get_data, "setup_cli_logging", _fake_logging)
    monkeypatch.setattr(get_data, "_warm_parent_catalog", lambda *_, **__: None)
    monkeypatch.setattr(get_data, "_ensure_testitem_pubchem_enabled", lambda *_, **__: None)
    monkeypatch.setattr(get_data, "_ensure_pubchem_enabled", lambda *_, **__: None)
    monkeypatch.setattr(get_data, "ensure_dirs", lambda *_, **__: None)

    class _ConfigStub:
        def __init__(self, base_dir: Path) -> None:
            self.base_dir = base_dir

        def model_copy(self, *, deep: bool = False) -> "_ConfigStub":
            return _ConfigStub(self.base_dir)

    monkeypatch.setattr(
        get_data,
        "_load_pipeline_config",
        lambda cfg, path: _ConfigStub(cfg.base_path),
    )

    def _fake_run_pipeline(cfg: Any, *, steps: Sequence[PipelineStep] | None = None) -> int:
        assert cfg.limit == 10
        effective_steps = tuple(steps or ())
        for step in effective_steps:
            spec = next(item for item in specs if item.name == step.name)
            payload = _HTTP_PAYLOADS[spec.endpoint]
            frame = pd.DataFrame(payload, columns=list(spec.columns)).astype("string")
            correlation = pd.DataFrame(
                {"metric": ["count"], "value": [str(len(frame))]}, dtype="string"
            )
            quality = pd.DataFrame(
                {
                    "column": frame.columns,
                    "non_null": frame.notna().sum().astype("int64").astype(str),
                }
            ).astype("string")
            save_standard_outputs(
                frame,
                correlation,
                quality,
                table_name=step.output_stem,
                date_tag=cfg.date_prefix,
                output_dir=cfg.output_dir,
            )
        return 0

    monkeypatch.setattr(get_data, "run_pipeline", _fake_run_pipeline)

    exit_code = get_data.main(
        [
            "--base-path",
            str(base_path),
            "--input-dir",
            "input",
            "--output-dir",
            "output",
            "--config",
            str(config_path),
            "--date",
            "20240115",
            "--limit",
            "10",
        ]
    )

    assert exit_code == 0

    produced = sorted(output_dir.glob("*.csv"))
    assert len(produced) == 15

    dataset_path = output_dir / "output.smoke_documents_20240115.csv"
    frame = pd.read_csv(dataset_path, dtype="string")
    assert list(frame.columns) == ["identifier", "title", "topic"]
    assert pd.isna(frame.loc[1, "topic"])
