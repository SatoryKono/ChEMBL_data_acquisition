from __future__ import annotations

from pathlib import Path

import pytest

from library.config import Config
from library.pipelines.assay import AssayPipelineOptions
from library.pipelines.assay import run_pipeline as run_assay_pipeline
from library.pipelines.common import PipelineRunResult
from library.pipelines.document import (
    DocumentPipelineOptions,
)
from library.pipelines.document import (
    run_pipeline as run_document_pipeline,
)
from library.pipelines.document import service as document_service
from library.pipelines.document.service import DocumentPipeline
from library.pipelines.target import TargetPipelineOptions
from library.pipelines.target import run_pipeline as run_target_pipeline
from scripts import get_assay_data, get_document_data, get_target_data


@pytest.fixture()
def sample_csv(tmp_path: Path) -> Path:
    path = tmp_path / "input.csv"
    path.write_text("id\nCHEMBL1\n", encoding="utf-8")
    return path


def test_run_document_service__invokes_mode_handler(
    cfg: Config, tmp_path: Path, sample_csv: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / ".output.documents_20240101.csv.tmp"
    cfg.io.output_dir = tmp_path
    options = DocumentPipelineOptions(
        input_csv=sample_csv,
        output_csv=output_csv,
        mode="chembl",
        offset=5,
        date_prefix="20240101",
        output_stem="documents",
    )

    captured: dict[str, object] = {}

    def _fake_run(
        cfg_arg: Config,
        args_arg,
        *,
        pipeline: DocumentPipeline | None = None,
    ) -> int:
        captured["cfg"] = cfg_arg
        captured["args"] = args_arg
        captured["pipeline"] = pipeline
        Path(args_arg.final_out).write_text("id\nCHEMBL1\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(get_document_data, "run_chembl", _fake_run)
    monkeypatch.setitem(get_document_data.MODE_HANDLERS, "chembl", _fake_run)

    monkeypatch.setattr(document_service, "_MODE_HANDLERS_CACHE", None)

    result = get_document_data.run_document_service(cfg, options)

    assert isinstance(result, PipelineRunResult)
    assert result.exit_code == 0
    assert result.output_path == output_csv
    assert result.executed is True
    assert result.written is True
    assert output_csv.exists()

    cfg_copy = captured["cfg"]
    assert isinstance(cfg_copy, Config)
    assert cfg_copy is not cfg
    assert (
        cfg_copy.sources.chembl.pipelines.document.chembl.offset
        == options.offset
    )

    args = captured["args"]
    assert Path(args.input_csv) == sample_csv
    assert Path(args.final_out) == tmp_path / "output.documents_20240101.csv"
    assert getattr(args, "_standard_date_tag") == "20240101"
    assert args.command == "chembl"

    pipeline = captured["pipeline"]
    assert isinstance(pipeline, DocumentPipeline)


def test_run_document_service__skip_existing(
    cfg: Config, tmp_path: Path, sample_csv: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / ".output.documents_20240101.csv.tmp"
    output_csv.write_text("existing", encoding="utf-8")
    options = DocumentPipelineOptions(
        input_csv=sample_csv,
        output_csv=output_csv,
        mode="chembl",
        skip_existing=True,
        date_prefix="20240101",
        output_stem="documents",
    )

    def _fail(*args, **kwargs):
        raise AssertionError("should not run")

    monkeypatch.setattr(get_document_data, "run_chembl", _fail)
    monkeypatch.setitem(get_document_data.MODE_HANDLERS, "chembl", _fail)

    monkeypatch.setattr(document_service, "_MODE_HANDLERS_CACHE", None)

    result = get_document_data.run_document_service(cfg, options)

    assert result.exit_code == 0
    assert result.executed is False
    assert result.reason == "skip_existing"
    assert result.written is False


def test_run_document_service__copies_canonical_output(
    cfg: Config, tmp_path: Path, sample_csv: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    working_output = tmp_path / ".output.documents_20240101.csv.tmp"
    cfg.io.output_dir = tmp_path
    options = DocumentPipelineOptions(
        input_csv=sample_csv,
        output_csv=working_output,
        mode="all",
        date_prefix="20240101",
        output_stem="documents",
    )

    def _fake_run(
        cfg_arg: Config,
        args_arg,
        *,
        pipeline: DocumentPipeline | None = None,
    ) -> int:
        Path(args_arg.final_out).write_text("id\nCHEMBL1\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(get_document_data, "run_all", _fake_run)
    monkeypatch.setitem(get_document_data.MODE_HANDLERS, "all", _fake_run)
    monkeypatch.setattr(document_service, "_MODE_HANDLERS_CACHE", None)

    result = get_document_data.run_document_service(cfg, options)

    assert result.exit_code == 0
    assert result.output_path == working_output
    assert working_output.exists()
    assert working_output.read_text(encoding="utf-8").splitlines()[1] == "CHEMBL1"


def test_run_document_service__copies_cli_output_when_canonical_missing(
    cfg: Config, tmp_path: Path, sample_csv: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    working_output = tmp_path / ".output.documents_20240101.csv.tmp"
    cfg.io.output_dir = tmp_path
    options = DocumentPipelineOptions(
        input_csv=sample_csv,
        output_csv=working_output,
        mode="chembl",
        date_prefix="20240101",
        output_stem="documents",
    )

    def _fake_run(
        cfg_arg: Config,
        args_arg,
        *,
        pipeline: DocumentPipeline | None = None,
    ) -> int:
        final_out = Path(args_arg.final_out)
        final_out.write_text("id\nCHEMBL2\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(get_document_data, "run_chembl", _fake_run)
    monkeypatch.setitem(get_document_data.MODE_HANDLERS, "chembl", _fake_run)
    monkeypatch.setattr(document_service, "_MODE_HANDLERS_CACHE", None)

    result = get_document_data.run_document_service(cfg, options)

    assert result.exit_code == 0
    assert result.output_path == working_output
    assert working_output.exists()
    assert working_output.read_text(encoding="utf-8").splitlines()[1] == "CHEMBL2"


def test_document_pipeline_run__delegates_to_service(
    cfg: Config, tmp_path: Path, sample_csv: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / "output.csv"
    options = DocumentPipelineOptions(
        input_csv=sample_csv,
        output_csv=output_csv,
        mode="chembl",
    )

    sentinel = PipelineRunResult(exit_code=0, output_path=output_csv, written=True)

    def _fake_service(config: Config, opts: DocumentPipelineOptions) -> PipelineRunResult:
        assert config is cfg
        assert opts is options
        return sentinel

    monkeypatch.setattr(document_service, "_MODE_HANDLERS_CACHE", None)
    monkeypatch.setattr(document_service, "run_document_service", _fake_service)

    result = run_document_pipeline(cfg, options)
    assert result is sentinel


def test_run_target_service__invokes_command_handler(
    cfg: Config, tmp_path: Path, sample_csv: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / "targets.csv"
    options = TargetPipelineOptions(
        input_csv=sample_csv,
        output_csv=output_csv,
        command="chembl",
        limit=10,
        offset=3,
        date="20250115",
        output_stem="custom_targets",
    )

    captured: dict[str, object] = {}

    def _fake_run(cfg_arg: Config, args_arg) -> int:
        captured["cfg"] = cfg_arg
        captured["args"] = args_arg
        return 0

    monkeypatch.setattr(get_target_data, "run_chembl", _fake_run)

    result = get_target_data.run_target_service(cfg, options)

    assert result.exit_code == 0
    assert result.output_path == output_csv
    assert result.executed is True
    assert result.written is True

    cfg_copy = captured["cfg"]
    assert isinstance(cfg_copy, Config)
    assert cfg_copy is not cfg
    assert (
        cfg_copy.sources.chembl.pipelines.target.chembl.limit
        == options.limit
    )
    assert (
        cfg_copy.sources.chembl.pipelines.target.chembl.offset
        == options.offset
    )

    args = captured["args"]
    assert Path(args.input_csv) == sample_csv
    assert Path(args.final_out) == output_csv
    assert args.command == "chembl"
    assert getattr(args, "_auto_output_stem", None) == options.output_stem
    assert getattr(args, "date", None) == options.date


def test_run_target_service__skip_existing(
    cfg: Config, tmp_path: Path, sample_csv: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / "targets.csv"
    output_csv.write_text("existing", encoding="utf-8")
    options = TargetPipelineOptions(
        input_csv=sample_csv,
        output_csv=output_csv,
        command="chembl",
        skip_existing=True,
    )

    monkeypatch.setattr(
        get_target_data,
        "run_chembl",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("should not run")),
    )

    result = get_target_data.run_target_service(cfg, options)

    assert result.exit_code == 0
    assert result.executed is False
    assert result.reason == "skip_existing"
    assert result.written is False


def test_target_pipeline_run__delegates_to_service(
    cfg: Config, tmp_path: Path, sample_csv: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / "targets.csv"
    options = TargetPipelineOptions(
        input_csv=sample_csv,
        output_csv=output_csv,
        command="chembl",
    )

    sentinel = PipelineRunResult(exit_code=0, output_path=output_csv, written=True)

    def _fake_service(config: Config, opts: TargetPipelineOptions) -> PipelineRunResult:
        assert config is cfg
        assert opts is options
        return sentinel

    monkeypatch.setattr(get_target_data, "run_target_service", _fake_service)

    result = run_target_pipeline(cfg, options)
    assert result is sentinel


def test_run_assay_service__invokes_cli_runner(
    cfg: Config, tmp_path: Path, sample_csv: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / "assays.csv"
    options = AssayPipelineOptions(
        input_csv=sample_csv,
        output_csv=output_csv,
        limit=4,
        offset=2,
        timeout=12.5,
        batch_size=7,
    )

    captured: dict[str, object] = {}

    def _fake_run(cfg_arg: Config, args_arg) -> int:
        captured["cfg"] = cfg_arg
        captured["args"] = args_arg
        return 0

    monkeypatch.setattr(get_assay_data, "run", _fake_run)

    result = get_assay_data.run_assay_service(cfg, options)

    assert result.exit_code == 0
    assert result.output_path == output_csv
    assert result.executed is True
    assert result.written is True

    cfg_copy = captured["cfg"]
    assert isinstance(cfg_copy, Config)
    assert cfg_copy is not cfg
    assert cfg_copy.sources.chembl.pipelines.assay.offset == options.offset
    assert cfg_copy.sources.chembl.pipelines.assay.limit == options.limit
    assert cfg_copy.sources.chembl.pipelines.assay.timeout == options.timeout
    assert cfg_copy.sources.chembl.pipelines.assay.batch_size == options.batch_size

    args = captured["args"]
    assert Path(args.input_csv) == sample_csv
    assert Path(args.final_out) == output_csv
    assert args.limit == options.limit
    assert args.offset == options.offset
    assert args.timeout == options.timeout
    assert args.batch_size == options.batch_size


def test_run_assay_service__skip_existing(
    cfg: Config, tmp_path: Path, sample_csv: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / "assays.csv"
    output_csv.write_text("existing", encoding="utf-8")
    options = AssayPipelineOptions(
        input_csv=sample_csv,
        output_csv=output_csv,
        skip_existing=True,
    )

    monkeypatch.setattr(
        get_assay_data,
        "run",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("should not run")),
    )

    result = get_assay_data.run_assay_service(cfg, options)

    assert result.exit_code == 0
    assert result.executed is False
    assert result.reason == "skip_existing"
    assert result.written is False


def test_assay_pipeline_run__delegates_to_service(
    cfg: Config, tmp_path: Path, sample_csv: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / "assays.csv"
    options = AssayPipelineOptions(
        input_csv=sample_csv,
        output_csv=output_csv,
    )

    sentinel = PipelineRunResult(exit_code=0, output_path=output_csv, written=True)

    def _fake_service(config: Config, opts: AssayPipelineOptions) -> PipelineRunResult:
        assert config is cfg
        assert opts is options
        return sentinel

    monkeypatch.setattr(get_assay_data, "run_assay_service", _fake_service)

    result = run_assay_pipeline(cfg, options)
    assert result is sentinel
