"""Integration tests for the CLI entrypoints in ``scripts/get_*`` modules."""

from __future__ import annotations

import argparse
from collections import Counter
from collections.abc import Callable, Iterable
from pathlib import Path

import pandas as pd
import pytest
import requests

from library import cli_utils
from library.cli import LoggerConfig
from library.cli import utils as cli_runtime_utils
from library.cli.base import PipelineCLIBase
from library.config import Config
from library.maintenance.legacy_outputs import SENTINEL_FILENAME
from library.pipelines.common import PipelineRunResult
from library.reporting.run_manifest import PipelineOutputReport
from library.utils.cli_tools import get_document_type
from scripts import (
    get_activity_data,
    get_assay_data,
    get_data,
    get_document_data,
    get_target_data,
    get_testitem_data,
)
from tests.helpers import ASSAY_ENRICHMENT_MIN_RATIO


@pytest.fixture()
def sample_csv(tmp_path: Path) -> Callable[[str], Path]:
    """Return a helper copying CSV fixtures into ``tmp_path``."""

    data_dir = Path(__file__).resolve().parents[1] / "resources" / "pipeline_inputs"

    def _copy(name: str) -> Path:
        source = data_dir / f"{name}.csv"
        target = tmp_path / f"{name}.csv"
        target.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
        return target

    return _copy


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


class _MemoryLogger:
    """Collect structured log events emitted by the CLI modules."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def _record(
        self,
        level: str,
        event: str,
        args: tuple[object, ...],
        kwargs: dict[str, object],
    ) -> None:
        payload = dict(kwargs)
        if args:
            payload["args"] = args
        self.events.append((level, event, payload))

    def info(self, event: str, *args: object, **kwargs: object) -> None:
        self._record("info", event, args, dict(kwargs))

    def warning(self, event: str, *args: object, **kwargs: object) -> None:
        self._record("warning", event, args, dict(kwargs))

    def error(self, event: str, *args: object, **kwargs: object) -> None:
        self._record("error", event, args, dict(kwargs))

    def debug(self, event: str, *args: object, **kwargs: object) -> None:
        self._record("debug", event, args, dict(kwargs))


def _patch_logger(monkeypatch: pytest.MonkeyPatch, module: object) -> _MemoryLogger:
    logger = _MemoryLogger()
    monkeypatch.setattr(module, "logger", logger)
    return logger


@pytest.mark.e2e
@pytest.mark.parametrize(
    "module",
    [
        get_data,
        get_activity_data,
        get_assay_data,
        get_document_data,
        get_target_data,
        get_testitem_data,
    ],
)
def test_cli_wrappers__monkeypatch_parent_catalog(
    monkeypatch: pytest.MonkeyPatch, module: object
) -> None:
    """Ensure monkeypatching private helpers affects both wrapper and underlying module."""

    sentinel = object()
    monkeypatch.setattr(module, "_warm_parent_catalog", sentinel, raising=False)

    assert module._warm_parent_catalog is sentinel
    assert module._MODULE._warm_parent_catalog is sentinel


@pytest.mark.parametrize(
    "module",
    [get_assay_data, get_activity_data, get_target_data],
)
def test_cli_wrappers_delegate_to_cli_instance(
    monkeypatch: pytest.MonkeyPatch, module: object
) -> None:
    """Ensure module-level helpers forward to :class:`PipelineCLIBase` instances."""

    assert isinstance(module._CLI, PipelineCLIBase)
    sentinel_parser = ("parser", "log")
    sentinel_exit = object()
    monkeypatch.setattr(module._CLI, "build_parser", lambda: sentinel_parser)
    monkeypatch.setattr(module._CLI, "main", lambda argv=None: sentinel_exit)

    assert module.build_parser() is sentinel_parser
    assert module.main([]) is sentinel_exit


class _DummyChemblClient:
    """Minimal context manager used to stub :class:`ChemblClient`."""

    def __init__(
        self, *args, **kwargs
    ) -> None:  # pragma: no cover - interface compatibility
        pass

    def __enter__(self) -> _DummyChemblClient:  # pragma: no cover - trivial helper
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:  # pragma: no cover - trivial helper
        return False

    def close(self) -> None:  # pragma: no cover - interface compatibility
        return None


@pytest.fixture()
def activity_resource_dir(snapshot_resource: Path) -> Path:
    return snapshot_resource / "activity_pipeline"


def _configure_activity_cfg(cfg: Config) -> None:
    cfg.activity.limit = None
    cfg.activity.offset = 0
    cfg.activity.dry_run = False
    cfg.activity.batch_size = 5
    cfg.activity.workers = 1
    cfg.system.doc_quality.enable = False
    cfg.activity_enrichment.action_type.enabled = False
    cfg.activity_enrichment.action_type.log_missing = False
    cfg.activity_enrichment.activity_properties.enabled = False


def _install_activity_writer(
    monkeypatch: pytest.MonkeyPatch,
) -> list[dict[str, object]]:
    written: list[dict[str, object]] = []

    def _writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        *,
        key_cols: Iterable[str] | None,
        col_order: Iterable[str] | None,
        chunksize: int,
        sort_chunksize: int | None = None,
        sep: str = ",",
        encoding: str = "utf-8",
        cfg=None,
        **_: object,
    ) -> Path:
        frames = [chunk.copy() for chunk in chunks]
        if frames:
            result = pd.concat(frames, ignore_index=True)
        else:
            result = pd.DataFrame(columns=list(col_order or []))
        if col_order:
            order = [str(col) for col in col_order]
            result = result.reindex(columns=order, fill_value=pd.NA)
        destination_path = Path(destination)
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(destination_path, index=False, sep=sep, encoding=encoding)

        quality_df = pd.DataFrame(
            [
                {
                    "metric": "rows_total",
                    "value": int(result.shape[0]),
                }
            ]
        )
        correlation_df = pd.DataFrame(
            [
                {
                    "column_a": "activity_id",
                    "column_b": "standard_value",
                    "correlation": 0.0,
                }
            ]
        )
        quality_path = destination_path.with_name(
            f"{destination_path.stem}_quality_report_table.csv"
        )
        correlation_path = destination_path.with_name(
            f"{destination_path.stem}_data_correlation_report_table.csv"
        )
        quality_df.to_csv(quality_path, index=False, sep=sep, encoding=encoding)
        correlation_df.to_csv(correlation_path, index=False, sep=sep, encoding=encoding)

        written.append(
            {
                "dataset_path": destination_path,
                "dataset": result,
                "quality_path": quality_path,
                "quality": quality_df,
                "correlation_path": correlation_path,
                "correlation": correlation_df,
                "sep": sep,
            }
        )
        return destination_path

    monkeypatch.setattr(get_activity_data, "write_csv_chunks_deterministic", _writer)
    return written


def _fake_activity_pipeline(
    params: dict[str, object],
    captured: dict[str, int],
    output_csv: Path,
) -> PipelineRunResult:
    fetch_config = params["fetch_config"]
    fetch_chunk = params["fetch_chunk"]
    metadata_hooks = params["metadata_hooks"]
    definition_kwargs = params["definition_kwargs"]
    writer_config = params["writer_config"]

    captured["workers"] = fetch_config.workers
    captured["chunk_size"] = fetch_config.chunk_size

    frames = [
        fetch_chunk(list(chunk_ids))
        for chunk_ids in fetch_config.chunker(fetch_config.ids, fetch_config.chunk_size)
    ]
    result = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    for hook in metadata_hooks:
        result = hook(result)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    writer = writer_config.writer
    writer_kwargs = dict(writer_config.kwargs)
    writer(
        [result],
        output_csv,
        list(result.columns),
        definition_kwargs.get("key_columns", ()),
        **writer_kwargs,
    )

    stats_callback = definition_kwargs.get("stats_callback")
    if stats_callback is not None:
        stats_callback(
            {"rows_total": len(result), "rows_kept": len(result), "rows_dropped": 0}
        )

    table_quality = definition_kwargs.get("table_quality")
    if callable(table_quality):
        table_quality(output_csv)

    return PipelineRunResult(exit_code=0, output_path=output_csv, written=True)


def _patch_activity_cli(monkeypatch: pytest.MonkeyPatch, cfg: Config) -> None:
    def fake_apply_config_overrides(
        args,
        parser,
        config_path,
        mapping=None,
        *,
        base_parser=None,
    ) -> Config:
        args._config_metadata = None
        final_out = getattr(args, "final_out", None)
        if final_out is not None:
            final_path = Path(final_out)
            args.final_out = final_path
            args.output_csv = final_path
        elif hasattr(args, "output_csv") and args.output_csv is not None:
            output_path = Path(args.output_csv)
            args.final_out = output_path
            args.output_csv = output_path
        cfg.activity.batch_size = getattr(args, "batch_size", cfg.activity.batch_size)
        cfg.activity.limit = getattr(args, "limit", cfg.activity.limit)
        cfg.activity.offset = getattr(args, "offset", cfg.activity.offset)
        cfg.activity.dry_run = getattr(args, "dry_run", cfg.activity.dry_run)
        cfg.activity.timeout = getattr(args, "timeout", cfg.activity.timeout)
        workers = getattr(args, "workers", None)
        if workers is not None:
            cfg.activity.workers = workers
        return cfg

    monkeypatch.setattr(
        cli_utils, "apply_config_overrides", fake_apply_config_overrides
    )
    monkeypatch.setattr(cli_utils, "ensure_dirs", lambda _cfg: None)


@pytest.mark.e2e
def test_run_cli_command__legacy_cleanup_removes_extra_files(
    tmp_path: Path,
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Ensure CLI boilerplate removes legacy artefacts before running pipelines."""

    cfg.io.output_dir = tmp_path / "output"
    cfg.io.cache_dir = tmp_path / "cache"
    cfg.io.exist_ok = True
    cfg.io.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.io.cache_dir.mkdir(parents=True, exist_ok=True)

    legacy_meta = cfg.io.output_dir / "documents.meta.yaml"
    legacy_failures = cfg.io.output_dir / "documents_failure_cases.csv"
    legacy_meta.write_text("meta", encoding="utf-8")
    legacy_failures.write_text("id\n1\n", encoding="utf-8")

    sentinel_path = cfg.io.output_dir / SENTINEL_FILENAME
    assert not sentinel_path.exists()

    logger_stub = _MemoryLogger()
    monkeypatch.setattr(
        cli_runtime_utils,
        "configure_logger",
        lambda _cfg: logger_stub,
    )

    def _fake_apply_config_overrides(
        args: argparse.Namespace,
        parser: argparse.ArgumentParser,
        config_path: Path | str,
        mapping: dict[str, str],
        *,
        base_parser: argparse.ArgumentParser | None = None,
    ) -> Config:
        args._config_metadata = None
        return cfg

    monkeypatch.setattr(
        cli_runtime_utils,
        "apply_config_overrides",
        _fake_apply_config_overrides,
    )
    monkeypatch.setattr(cli_runtime_utils, "prepare_io_paths", lambda _args: None)

    output_csv = cfg.io.output_dir / "documents.csv"

    run_calls: list[tuple[Config, argparse.Namespace]] = []

    def _run(cfg_obj: Config, args_ns: argparse.Namespace) -> int:
        run_calls.append((cfg_obj, args_ns))
        output_csv.write_text("document_id\nDOC1\n", encoding="utf-8")
        return 0

    config_path = tmp_path / "config.yaml"
    config_path.write_text("{}", encoding="utf-8")

    parser = argparse.ArgumentParser(prog="get-document-data")
    parser.add_argument("--config", default=str(config_path))

    args = argparse.Namespace(
        config=str(config_path),
        log_level="INFO",
        verbose=False,
        run_id=None,
        print_config=False,
        postprocess=False,
        base_path=None,
        input_dir=None,
        output_dir=cfg.io.output_dir,
        cache_dir=cfg.io.cache_dir,
        input_csv=None,
        output_csv=output_csv,
        final_out=output_csv,
        invocation=("get-document-data",),
    )

    log_cfg = LoggerConfig(level="INFO", run_id="test-run")

    exit_code = cli_runtime_utils.run_cli_command(
        args=args,
        parser=parser,
        log_cfg=log_cfg,
        mapping={},
        run=_run,
    )

    assert exit_code == 0
    assert run_calls, "expected pipeline run to be invoked"
    assert output_csv.exists()
    assert not legacy_meta.exists()
    assert not legacy_failures.exists()
    assert sentinel_path.exists()

    events = [(level, event) for level, event, _ in logger_stub.events]
    assert ("warning", "legacy_outputs_retention_notice") in events


@pytest.mark.e2e
def test_get_activity_cli__default_date_prefix_applied(
    tmp_path: Path,
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Ensure CLI derives filenames and metadata from configured date prefix."""

    _configure_activity_cfg(cfg)
    cfg.io.output_dir = tmp_path / "out"
    cfg.io.cache_dir = tmp_path / "cache"
    cfg.io.default_date_prefix = "19991231"
    cfg.io.output_stamp_mode = "omit"

    input_csv = tmp_path / "activities.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")

    written = _install_activity_writer(monkeypatch)

    original_write_meta = get_activity_data.write_meta_yaml
    captured_meta: dict[str, object] = {}

    def _record_meta(
        path: Path | str, *, cfg: Config | None = None, **kwargs: object
    ) -> Path:
        meta_path = original_write_meta(path, cfg=cfg, **kwargs)
        captured_meta["path"] = meta_path
        captured_meta["cfg"] = cfg
        return meta_path

    monkeypatch.setattr(get_activity_data, "write_meta_yaml", _record_meta)

    def _stub_activity_pipeline(**kwargs: object) -> PipelineRunResult:
        output_csv = Path(kwargs["output_path"])
        writer_cfg = kwargs["writer_config"]
        definition_kwargs = dict(kwargs.get("definition_kwargs", {}))
        frame = pd.DataFrame([{"activity_id": "ACT1"}])
        writer_cfg.writer(
            [frame],
            output_csv,
            key_cols=definition_kwargs.get("key_columns", ()),
            col_order=definition_kwargs.get("column_order", tuple(frame.columns)),
            **writer_cfg.kwargs,
        )
        get_activity_data.write_meta_yaml(
            output_csv,
            cfg=kwargs["cfg"],
            columns=list(frame.columns),
        )
        return PipelineRunResult(exit_code=0, output_path=output_csv, written=True)

    monkeypatch.setattr(
        "library.pipelines.activity.run.run_activity_pipeline",
        _stub_activity_pipeline,
    )

    _patch_activity_cli(monkeypatch, cfg)

    exit_code = get_activity_data.main(["--input", str(input_csv)])

    assert exit_code == 0
    assert len(written) == 1

    output_path = written[0][0]
    expected_name = f"output.{get_activity_data.DEFAULT_OUTPUT_STEM}_{cfg.io.default_date_prefix}.csv"
    assert output_path.name == expected_name
    assert output_path.parent == cfg.io.output_dir

    meta_path = captured_meta.get("path")
    assert isinstance(meta_path, Path)
    assert meta_path.name == f"{expected_name}.meta.yaml"
    assert meta_path.parent == output_path.parent
    assert captured_meta.get("cfg") is cfg


@pytest.mark.e2e
def test_get_document_type_main__writes_meta(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "docs.csv"
    frame = pd.DataFrame(
        [
            {
                "chembl_id": "CHEMBL1",
                "PubMed.PublicationType": "Review Article",
                "scholar.PublicationTypes": "",
                "OpenAlex.PublicationTypes": "",
            },
            {
                "chembl_id": "CHEMBL2",
                "PubMed.PublicationType": "Journal Article",
                "scholar.PublicationTypes": "Conference Paper",
                "OpenAlex.PublicationTypes": "",
            },
        ]
    )
    frame.to_csv(input_csv, index=False)
    output_csv = tmp_path / "output" / "document_types.csv"

    def fake_apply_config_overrides(
        args,
        parser,
        config_path,
        mapping=None,
        *,
        base_parser=None,
    ) -> Config:
        args._config_metadata = None
        if getattr(args, "output_csv", None) is not None:
            output_path = Path(args.output_csv)
            args.output_csv = output_path
            args.final_out = output_path
        cfg.io.output_dir = tmp_path / "io"
        cfg.io.cache_dir = tmp_path / "cache"
        return cfg

    monkeypatch.setattr(
        get_document_type.cli, "apply_config_overrides", fake_apply_config_overrides
    )

    captured: dict[str, object] = {}
    finalise_payload: dict[str, object] = {}

    def _stub_write_csv(
        frame: pd.DataFrame,
        path: Path,
        *,
        cfg: Config,
        key_cols: Iterable[str] | None = None,
        sep: str | None = None,
        encoding: str | None = None,
        **_: object,
    ) -> Path:
        dataset = frame.copy()
        dataset_path = Path(path)
        dataset_path.parent.mkdir(parents=True, exist_ok=True)
        dataset.to_csv(
            dataset_path,
            index=False,
            sep=sep or ",",
            encoding=encoding or "utf-8",
        )

        quality_df = pd.DataFrame(
            [
                {
                    "metric": "rows_total",
                    "value": int(dataset.shape[0]),
                }
            ]
        )
        correlation_df = pd.DataFrame(
            [
                {
                    "column_a": "chembl_id",
                    "column_b": "class_label",
                    "correlation": 1.0,
                }
            ]
        )
        quality_path = dataset_path.with_name(
            f"{dataset_path.stem}_quality_report_table.csv"
        )
        correlation_path = dataset_path.with_name(
            f"{dataset_path.stem}_data_correlation_report_table.csv"
        )
        quality_df.to_csv(
            quality_path,
            index=False,
            sep=sep or ",",
            encoding=encoding or "utf-8",
        )
        correlation_df.to_csv(
            correlation_path,
            index=False,
            sep=sep or ",",
            encoding=encoding or "utf-8",
        )

        captured["cfg"] = cfg
        captured["key_cols"] = list(key_cols or [])
        captured["dataset"] = dataset
        captured["quality_df"] = quality_df
        captured["correlation_df"] = correlation_df
        captured["paths"] = {
            "dataset": dataset_path,
            "quality": quality_path,
            "correlation": correlation_path,
        }
        return dataset_path

    monkeypatch.setattr(get_document_type.io, "write_csv", _stub_write_csv)
    monkeypatch.setattr(get_document_type, "write_meta_yaml", lambda *_, **__: None)

    def _fake_finalise(**kwargs):
        finalise_payload["kwargs"] = kwargs
        csv_path = Path(kwargs["csv_path"])
        stats = {
            "rows_total": kwargs.get("rows_total", 0),
            "rows_kept": kwargs.get("rows_kept", 0),
            "rows_dropped": kwargs.get("rows_total", 0) - kwargs.get("rows_kept", 0),
            "output_sha256": "stub",
        }
        return PipelineOutputReport(
            csv_path=csv_path,
            stats=stats,
            meta_path=csv_path.with_name(csv_path.name + ".meta.yaml"),
            meta_sha256="meta",
        )

    monkeypatch.setattr("library.cli.utils.finalise_csv_output", _fake_finalise)

    invocation = ["--input", str(input_csv), "--final-out", str(output_csv)]
    exit_code = get_document_type.main(invocation)

    assert exit_code == 0
    assert captured["cfg"] is cfg
    assert captured["key_cols"] == ["chembl_id"]

    assert output_csv.exists()
    result = pd.read_csv(output_csv)
    assert "class_label" in result.columns
    assert captured["paths"]["dataset"] == output_csv
    quality_path = captured["paths"]["quality"]
    correlation_path = captured["paths"]["correlation"]
    assert quality_path.exists()
    assert correlation_path.exists()

    produced_csvs = {path.name for path in output_csv.parent.glob("*.csv")}
    assert produced_csvs == {
        output_csv.name,
        quality_path.name,
        correlation_path.name,
    }
    assert not list(output_csv.parent.glob("*.meta.yaml"))

    finalise_kwargs = finalise_payload.get("kwargs") or {}
    assert finalise_kwargs.get("schema") == "document_type"
    assert finalise_kwargs.get("rows_total") == len(frame)
    assert finalise_kwargs.get("rows_kept") == len(frame)
    assert finalise_kwargs.get("command")
    assert finalise_kwargs.get("invocation") == tuple(invocation)


@pytest.mark.e2e
def test_get_testitem_run_success(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("testitem")
    output_csv = tmp_path / "out" / "testitem.csv"
    logger_stub = _patch_logger(monkeypatch, get_testitem_data)

    def _stub_pipeline(
        config: Config, options: get_testitem_data.TestitemPipelineOptions
    ) -> int:
        assert config is cfg
        assert options.input_csv == Path(input_csv)
        assert options.output_csv == output_csv
        frame = pd.read_csv(options.input_csv)
        frame["preferred_name"] = (
            frame["preferred_name"].fillna("").astype("string").str.strip()
        )
        missing = frame["preferred_name"] == ""
        if int(missing.sum()):
            get_testitem_data.logger.warning(
                "testitem_missing_name", count=int(missing.sum())
            )
        frame["normalized_name"] = (
            frame["preferred_name"].replace("", pd.NA).str.lower()
        )
        frame["is_named"] = (~missing).astype("boolean")
        _ensure_parent(output_csv)
        frame.to_csv(output_csv, index=False)
        return 0

    monkeypatch.setattr(get_testitem_data, "run_testitem_pipeline", _stub_pipeline)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
    )

    rc = get_testitem_data.run(cfg, args)

    assert rc == 0
    artifacts = getattr(args, "_testitem_artifacts", None)
    assert artifacts is not None, "expected standard output artifacts"
    dataset_path = artifacts.dataset
    result = pd.read_csv(dataset_path)
    assert list(result.columns) == [
        "molecule_chembl_id",
        "preferred_name",
        "normalized_name",
        "is_named",
    ]
    assert result.loc[0, "normalized_name"] == "compound 1"
    assert pd.isna(result.loc[1, "normalized_name"])
    assert not result.loc[1, "is_named"]
    meta_path = dataset_path.with_suffix(".meta.yaml")
    assert meta_path.exists()
    events = [event for _, event, _ in logger_stub.events]
    assert "testitem_pipeline_done" in events
    assert "testitem_missing_name" in events


@pytest.mark.e2e
def test_get_testitem_run_failure_logs(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("testitem")
    output_csv = tmp_path / "out" / "testitem.csv"
    logger_stub = _patch_logger(monkeypatch, get_testitem_data)

    def _failing_pipeline(
        config: Config, options: get_testitem_data.TestitemPipelineOptions
    ) -> int:
        get_testitem_data.logger.error(
            "testitem_pipeline_failed", output=str(options.output_csv), exit_code=2
        )
        return 2

    monkeypatch.setattr(get_testitem_data, "run_testitem_pipeline", _failing_pipeline)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
    )

    rc = get_testitem_data.run(cfg, args)

    assert rc == 2
    events = [event for _, event, _ in logger_stub.events]
    assert "testitem_pipeline_failed" in events


@pytest.mark.e2e
def test_get_testitem_run_skip_existing(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("testitem")
    output_csv = tmp_path / "out" / "testitem.csv"
    _ensure_parent(output_csv)
    output_csv.write_text("placeholder", encoding="utf-8")

    call_counter = Counter()
    logger_stub = _patch_logger(monkeypatch, get_testitem_data)

    def _unexpected_call(*_: object, **__: object) -> int:
        call_counter["called"] += 1
        return 0

    monkeypatch.setattr(get_testitem_data, "run_testitem_pipeline", _unexpected_call)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=True,
        force=False,
    )

    rc = get_testitem_data.run(cfg, args)

    assert rc == 0
    assert call_counter["called"] == 0
    events = [event for _, event, _ in logger_stub.events]
    assert "pipeline_skip_existing" in events


@pytest.mark.e2e
def test_get_activity_cli__retry_and_idempotent(
    tmp_path: Path,
    activity_resource_dir: Path,
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _configure_activity_cfg(cfg)
    input_csv = tmp_path / "activities_input.csv"
    input_csv.write_text(
        (activity_resource_dir / "ids_happy.csv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    output_csv = tmp_path / "out" / "activities.csv"
    chunk_df = pd.read_csv(activity_resource_dir / "chunk_happy.csv")

    attempts = {"count": 0}

    def _fake_get_activities(chunk_ids, **_kwargs):
        attempts["count"] += 1
        if attempts["count"] % 2 == 1:
            raise requests.RequestException("temporary failure")
        identifiers = [str(item) for item in chunk_ids]
        mask = chunk_df["activity_id"].astype(str).isin(identifiers)
        return chunk_df.loc[mask].reset_index(drop=True)

    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        _DummyChemblClient,
    )
    monkeypatch.setattr(get_activity_data.cl, "get_activities", _fake_get_activities)
    monkeypatch.setattr(get_activity_data, "sleep", lambda *_args, **_kwargs: None)

    written = _install_activity_writer(monkeypatch)
    logger_stub = _patch_logger(monkeypatch, get_activity_data)
    _patch_activity_cli(monkeypatch, cfg)

    args = ["--input", str(input_csv), "--final-out", str(output_csv)]

    first_exit = get_activity_data.main(args)
    assert first_exit == 0
    assert output_csv.exists()
    first_content = output_csv.read_text(encoding="utf-8")
    events = [event for _, event, _ in logger_stub.events]
    assert "activity_fetch_retry" in events
    assert "activity_pipeline_done" in events
    assert not any(event == "activity_pipeline_failed" for event in events)
    assert attempts["count"] >= 2
    assert len(written) == 1
    assert list(written[0]["dataset"]["activity_id"]) == ["ACT1", "ACT2", "ACT3"]

    logger_stub.events.clear()
    written.clear()
    attempts["count"] = 0

    second_exit = get_activity_data.main(args)
    assert second_exit == 0
    second_content = output_csv.read_text(encoding="utf-8")
    assert second_content == first_content
    done_events = [event for _, event, _ in logger_stub.events]
    assert done_events.count("activity_pipeline_done") >= 1
    assert not any(event == "activity_pipeline_failed" for event in done_events)


@pytest.mark.e2e
def test_get_activity_cli__timeout_split_recovers(
    tmp_path: Path,
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _configure_activity_cfg(cfg)
    cfg.activity.batch_size = 4
    cfg.retry.max_attempts = 2
    input_csv = tmp_path / "activities.csv"
    input_csv.write_text(
        "activity_id\nACT1\nACT2\nACT3\nACT4\n",
        encoding="utf-8",
    )
    output_csv = tmp_path / "out" / "activities.csv"

    activity_rows = {
        "ACT1": {
            "activity_id": "ACT1",
            "molecule_chembl_id": "CHEMBL1",
            "assay_chembl_id": "ASSAY1",
            "standard_value": 1.0,
            "standard_units": "nM",
            "standard_type": "IC50",
            "relation": "=",
            "molecule_pref_name": "Compound 1",
        },
        "ACT2": {
            "activity_id": "ACT2",
            "molecule_chembl_id": "CHEMBL2",
            "assay_chembl_id": "ASSAY2",
            "standard_value": 2.0,
            "standard_units": "nM",
            "standard_type": "IC50",
            "relation": "=",
            "molecule_pref_name": "Compound 2",
        },
        "ACT3": {
            "activity_id": "ACT3",
            "molecule_chembl_id": "CHEMBL3",
            "assay_chembl_id": "ASSAY3",
            "standard_value": 3.0,
            "standard_units": "nM",
            "standard_type": "IC50",
            "relation": "=",
            "molecule_pref_name": "Compound 3",
        },
        "ACT4": {
            "activity_id": "ACT4",
            "molecule_chembl_id": "CHEMBL4",
            "assay_chembl_id": "ASSAY4",
            "standard_value": 4.0,
            "standard_units": "nM",
            "standard_type": "IC50",
            "relation": "=",
            "molecule_pref_name": "Compound 4",
        },
    }

    call_log: list[tuple[str, ...]] = []

    def _fake_get_activities(chunk_ids: Iterable[str], **_: object) -> pd.DataFrame:
        identifiers = [str(item) for item in chunk_ids]
        call_log.append(tuple(identifiers))
        if len(identifiers) > 1:
            raise requests.ReadTimeout("simulated timeout")
        row = activity_rows[identifiers[0]]
        return pd.DataFrame([row])

    monkeypatch.setattr(get_activity_data.cl, "get_activities", _fake_get_activities)
    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        _DummyChemblClient,
    )
    monkeypatch.setattr(get_activity_data, "sleep", lambda *_args, **_kwargs: None)

    written = _install_activity_writer(monkeypatch)
    logger_stub = _patch_logger(monkeypatch, get_activity_data)
    _patch_activity_cli(monkeypatch, cfg)

    args = [
        "--input",
        str(input_csv),
        "--final-out",
        str(output_csv),
        "--batch-size",
        "4",
    ]

    exit_code = get_activity_data.main(args)

    assert exit_code == 0
    assert output_csv.exists()
    result = pd.read_csv(output_csv)
    assert list(result["activity_id"]) == ["ACT1", "ACT2", "ACT3", "ACT4"]
    assert len(written) == 1

    events = [event for _, event, _ in logger_stub.events]
    assert "activity_fetch_split" in events
    assert "activity_pipeline_done" in events
    assert not any(event == "activity_pipeline_failed" for event in events)

    assert call_log[0] == ("ACT1", "ACT2", "ACT3", "ACT4")
    assert any(len(ids) == 1 for ids in call_log)


@pytest.mark.e2e
def test_get_activity_cli__workers_and_offset(
    tmp_path: Path,
    activity_resource_dir: Path,
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _configure_activity_cfg(cfg)
    input_csv = tmp_path / "activities_input.csv"
    input_csv.write_text(
        (activity_resource_dir / "ids_happy.csv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    output_csv = tmp_path / "out" / "activities.csv"
    chunk_df = pd.read_csv(activity_resource_dir / "chunk_happy.csv")

    def _fake_get_activities(chunk_ids, **_kwargs):
        identifiers = [str(item) for item in chunk_ids]
        mask = chunk_df["activity_id"].astype(str).isin(identifiers)
        return chunk_df.loc[mask].reset_index(drop=True)

    captured: dict[str, int] = {}

    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        _DummyChemblClient,
    )
    monkeypatch.setattr(get_activity_data.cl, "get_activities", _fake_get_activities)
    monkeypatch.setattr(
        "library.pipelines.activity.run.run_activity_pipeline",
        lambda **kwargs: _fake_activity_pipeline(kwargs, captured, output_csv),
    )

    written = _install_activity_writer(monkeypatch)
    logger_stub = _patch_logger(monkeypatch, get_activity_data)
    _patch_activity_cli(monkeypatch, cfg)

    args = [
        "--input",
        str(input_csv),
        "--final-out",
        str(output_csv),
        "--workers",
        "2",
        "--offset",
        "1",
        "--batch-size",
        "2",
    ]

    exit_code = get_activity_data.main(args)

    assert exit_code == 0
    assert captured["workers"] == 2
    assert output_csv.exists()
    assert len(written) == 1
    written_df = written[0]["dataset"]
    assert written_df["activity_id"].tolist() == ["ACT2", "ACT3"]
    events = [event for _, event, _ in logger_stub.events]
    assert "process_offset" in events
    assert "activity_pipeline_done" in events


@pytest.mark.e2e
def test_get_activity_cli__chembl_identifier_backfill_ratio(
    tmp_path: Path,
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _configure_activity_cfg(cfg)

    total_rows = 200
    activity_ids = [f"ACT{i:05d}" for i in range(total_rows)]
    chunk_df = pd.DataFrame(
        {
            "activity_id": pd.Series(activity_ids, dtype="string"),
            "activity_chembl_id": pd.Series([pd.NA] * total_rows, dtype="string"),
            "molecule_chembl_id": pd.Series(
                [f"CHEMBL{i % 7:05d}" for i in range(total_rows)], dtype="string"
            ),
            "assay_chembl_id": pd.Series(
                [f"ASSAY{i % 11:05d}" for i in range(total_rows)], dtype="string"
            ),
            "standard_value": pd.Series([float(i % 50 + 1) for i in range(total_rows)]),
            "standard_units": pd.Series(["nM"] * total_rows, dtype="string"),
            "standard_type": pd.Series(["IC50"] * total_rows, dtype="string"),
            "relation": pd.Series(["="] * total_rows, dtype="string"),
        }
    )

    input_csv = tmp_path / "activities_input.csv"
    input_csv.write_text(
        "activity_id\n" + "\n".join(activity_ids),
        encoding="utf-8",
    )

    def _fake_get_activities(chunk_ids, **_kwargs):
        identifiers = [str(item) for item in chunk_ids]
        mask = chunk_df["activity_id"].astype(str).isin(identifiers)
        return chunk_df.loc[mask].reset_index(drop=True)

    output_csv = tmp_path / "out" / "activities.csv"

    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        _DummyChemblClient,
    )
    monkeypatch.setattr(get_activity_data.cl, "get_activities", _fake_get_activities)

    written = _install_activity_writer(monkeypatch)
    logger_stub = _patch_logger(monkeypatch, get_activity_data)
    _patch_activity_cli(monkeypatch, cfg)

    args = ["--input", str(input_csv), "--final-out", str(output_csv)]

    exit_code = get_activity_data.main(args)

    assert exit_code == 0
    assert output_csv.exists()
    assert written, "expected pipeline to emit output"

    written_df = written[0]["dataset"]
    assert "activity_id" in written_df.columns
    assert "activity_chembl_id" in written_df.columns

    id_series = written_df["activity_id"].astype("string")
    chembl_series = written_df["activity_chembl_id"].astype("string")
    id_present = id_series.notna() & id_series.str.strip().ne("")
    assert int(id_present.sum()) == total_rows

    chembl_present = chembl_series.notna() & chembl_series.str.strip().ne("")
    filled_ratio = chembl_present[id_present].sum() / id_present.sum()
    assert filled_ratio >= 0.99

    # Spot-check fallback alignment for a few rows.
    assert (
        chembl_series[id_present].iloc[:5].tolist()
        == id_series[id_present].iloc[:5].tolist()
    )

    events = [event for _, event, _ in logger_stub.events]
    assert "activity_pipeline_done" in events
    assert "activity_pipeline_failed" not in events


@pytest.mark.e2e
def test_get_activity_cli__non_csv_output_path(
    tmp_path: Path,
    activity_resource_dir: Path,
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _configure_activity_cfg(cfg)
    cfg.io.csv_sep = "\t"
    input_csv = tmp_path / "activities_input.csv"
    input_csv.write_text(
        (activity_resource_dir / "ids_happy.csv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    output_csv = tmp_path / "out" / "activities.tsv"
    chunk_df = pd.read_csv(activity_resource_dir / "chunk_happy.csv")

    def _fake_get_activities(chunk_ids, **_kwargs):
        identifiers = [str(item) for item in chunk_ids]
        mask = chunk_df["activity_id"].astype(str).isin(identifiers)
        return chunk_df.loc[mask].reset_index(drop=True)

    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        _DummyChemblClient,
    )
    monkeypatch.setattr(get_activity_data.cl, "get_activities", _fake_get_activities)

    written = _install_activity_writer(monkeypatch)
    logger_stub = _patch_logger(monkeypatch, get_activity_data)
    _patch_activity_cli(monkeypatch, cfg)

    exit_code = get_activity_data.main(
        ["--input", str(input_csv), "--final-out", str(output_csv)]
    )

    assert exit_code == 0
    assert output_csv.exists()
    content = output_csv.read_text(encoding="utf-8")
    assert "\t" in content.splitlines()[1]
    assert len(written) == 1
    assert written[0]["sep"] == "\t"
    events = [event for _, event, _ in logger_stub.events]
    assert "activity_pipeline_done" in events


@pytest.mark.e2e
def test_get_document_run_all_success(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("document")
    output_csv = tmp_path / "out" / "documents.csv"
    logger_stub = _patch_logger(monkeypatch, get_document_data)

    def _stub_all(config: Config, args: argparse.Namespace) -> int:
        assert config is cfg
        frame = pd.read_csv(args.input_csv)
        frame["title"] = frame["title"].astype("string").str.strip()
        frame["pubmed_id"] = frame["pubmed_id"].astype("string").str.strip()
        frame["has_pubmed"] = frame["pubmed_id"].replace("", pd.NA).notna()
        frame = frame.drop_duplicates(subset=["document_chembl_id"])
        frame = frame.sort_values("document_chembl_id").reset_index(drop=True)
        missing = (~frame["has_pubmed"]).sum()
        if int(missing):
            get_document_data.logger.warning(
                "document_missing_pubmed", count=int(missing)
            )
        output_path = Path(args.final_out)
        _ensure_parent(output_path)
        frame.to_csv(output_path, index=False)
        get_document_data.logger.info(
            "document_all_done", output=str(args.final_out), rows=len(frame)
        )
        return 0

    monkeypatch.setattr(get_document_data, "run_all", _stub_all)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
        command="all",
        func=_stub_all,
        timeout=None,
        chunk_size=5,
    )

    rc = get_document_data.run(cfg, args)

    assert rc == 0
    result = pd.read_csv(output_csv)
    assert list(result.columns) == [
        "document_chembl_id",
        "title",
        "pubmed_id",
        "has_pubmed",
    ]
    assert set(result["document_chembl_id"]) == {"CHEMBL123", "CHEMBL456"}
    assert not result.loc[
        result["document_chembl_id"] == "CHEMBL456", "has_pubmed"
    ].iloc[0]
    events = [event for _, event, _ in logger_stub.events]
    assert "document_all_done" in events
    assert "document_missing_pubmed" in events


@pytest.mark.e2e
def test_get_document_run_missing_handler(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("document")
    output_csv = tmp_path / "out" / "documents.csv"
    logger_stub = _patch_logger(monkeypatch, get_document_data)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
        command="all",
        func=None,
    )

    rc = get_document_data.run(cfg, args)

    assert rc == 1
    events = [event for _, event, _ in logger_stub.events]
    assert "missing_subcommand_handler" in events


@pytest.mark.e2e
def test_get_document_run_all_failure(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("document")
    output_csv = tmp_path / "out" / "documents.csv"
    logger_stub = _patch_logger(monkeypatch, get_document_data)

    def _failing_all(config: Config, args: argparse.Namespace) -> int:
        get_document_data.logger.error(
            "document_all_failed", output=str(args.final_out), exit_code=1
        )
        return 1

    monkeypatch.setattr(get_document_data, "run_all", _failing_all)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
        command="all",
        func=_failing_all,
    )

    rc = get_document_data.run(cfg, args)

    assert rc == 1
    events = [event for _, event, _ in logger_stub.events]
    assert "document_all_failed" in events


@pytest.mark.e2e
def test_get_target_run_all_success(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("target")
    final_out = tmp_path / "out" / "targets.csv"
    logger_stub = _patch_logger(monkeypatch, get_target_data)

    def _stub_all(config: Config, args: argparse.Namespace) -> int:
        frame = pd.read_csv(args.input_csv)
        frame["target_name"] = frame["target_name"].astype("string").str.strip()
        frame["organism"] = frame["organism"].astype("string").str.strip()
        frame["name_upper"] = frame["target_name"].str.upper()
        frame["has_organism"] = frame["organism"].replace("", pd.NA).notna()
        frame = frame.sort_values("target_chembl_id").reset_index(drop=True)
        _ensure_parent(args.final_out)
        frame.to_csv(args.final_out, index=False)
        get_target_data.logger.info(
            "target_all_done", output=str(args.final_out), rows=len(frame)
        )
        return 0

    monkeypatch.setattr(get_target_data, "run_all", _stub_all)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=final_out,
        raw_out=None,
        raw_format="csv",
        no_reindex_raw=False,
        skip_existing=False,
        force=False,
        command="all",
        func=_stub_all,
        chunk_size=10,
        batch_size=100,
        offset=0,
        id_cols=None,
    )

    rc = get_target_data.run(cfg, args)

    assert rc == 0
    result = pd.read_csv(final_out)
    assert list(result.columns) == [
        "target_chembl_id",
        "target_name",
        "organism",
        "name_upper",
        "has_organism",
    ]
    assert (result["name_upper"] == result["target_name"].str.upper()).all()
    events = [event for _, event, _ in logger_stub.events]
    assert "target_all_done" in events


@pytest.mark.e2e
def test_get_target_run_skip_existing(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("target")
    final_out = tmp_path / "out" / "targets.csv"
    _ensure_parent(final_out)
    final_out.write_text("existing", encoding="utf-8")
    logger_stub = _patch_logger(monkeypatch, get_target_data)

    def _unexpected_call(*_: object, **__: object) -> int:
        raise AssertionError("run_all should not be invoked when skipping")

    monkeypatch.setattr(get_target_data, "run_all", _unexpected_call)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=final_out,
        raw_out=None,
        raw_format="csv",
        no_reindex_raw=False,
        skip_existing=True,
        force=False,
        command="all",
        func=_unexpected_call,
        chunk_size=10,
        batch_size=100,
        offset=0,
        id_cols=None,
    )

    rc = get_target_data.run(cfg, args)

    assert rc == 0
    events = [event for _, event, _ in logger_stub.events]
    assert "pipeline_skip_existing" in events


@pytest.mark.e2e
def test_get_target_run_all_failure(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("target")
    final_out = tmp_path / "out" / "targets.csv"
    logger_stub = _patch_logger(monkeypatch, get_target_data)

    def _failing_all(config: Config, args: argparse.Namespace) -> int:
        get_target_data.logger.error(
            "pipeline_step_failed", step="all", output=str(args.final_out)
        )
        return 1

    monkeypatch.setattr(get_target_data, "run_all", _failing_all)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=final_out,
        raw_out=None,
        raw_format="csv",
        no_reindex_raw=False,
        skip_existing=False,
        force=False,
        command="all",
        func=_failing_all,
        chunk_size=10,
        batch_size=100,
        offset=0,
        id_cols=None,
    )

    rc = get_target_data.run(cfg, args)

    assert rc == 1
    events = [event for _, event, _ in logger_stub.events]
    assert "pipeline_step_failed" in events


@pytest.mark.e2e
def test_get_assay_run_success(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("assay")
    output_csv = tmp_path / "out" / "assays.csv"
    logger_stub = _patch_logger(monkeypatch, get_assay_data)
    data_dir = Path(__file__).resolve().parents[1] / "resources" / "pipeline_inputs"
    dictionary_path = data_dir / "assay_dictionary.csv"

    def _stub_run_chembl(config: Config, args: argparse.Namespace) -> int:
        frame = pd.read_csv(args.input_csv)
        dictionary = pd.read_csv(dictionary_path)
        dictionary["assay_chembl_id"] = dictionary["assay_chembl_id"].astype("string")
        enriched = frame.merge(dictionary, on="assay_chembl_id", how="left")
        enriched["description"] = enriched["description"].astype("string").str.strip()
        enriched["description_length"] = (
            enriched["description"].str.len().astype("Int64")
        )
        enriched["year"] = pd.to_numeric(enriched["year"], errors="coerce").astype(
            "Int64"
        )
        enriched = enriched.drop_duplicates(subset=["assay_chembl_id"])
        enriched = enriched.sort_values("assay_chembl_id").reset_index(drop=True)
        quality_columns = ["assay_strain", "assay_group", "year", "accession"]
        completeness = 1.0 - enriched[quality_columns].isna().mean()
        if float(completeness.min()) < ASSAY_ENRICHMENT_MIN_RATIO:
            raise AssertionError(
                "assay enrichment below threshold "
                f"(threshold={ASSAY_ENRICHMENT_MIN_RATIO}, completeness={completeness.to_dict()})"
            )
        output_path = Path(args.final_out)
        _ensure_parent(output_path)
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
        enriched.to_csv(output_path, index=False, columns=columns)
        get_assay_data.logger.info(
            "assay_pipeline_done",
            output_postprocessed=str(args.final_out),
            processed=len(enriched),
        )
        return 0

    monkeypatch.setattr(get_assay_data, "run_chembl", _stub_run_chembl)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
    )

    rc = get_assay_data.run(cfg, args)

    assert rc == 0
    result = pd.read_csv(output_csv)
    assert list(result.columns) == [
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
    assert (result["description_length"] == result["description"].str.len()).all()
    quality_columns = ["assay_strain", "assay_group", "year", "accession"]
    completeness = 1.0 - result[quality_columns].isna().mean()
    assert (completeness >= ASSAY_ENRICHMENT_MIN_RATIO).all(), completeness
    events = [event for _, event, _ in logger_stub.events]
    assert "assay_pipeline_done" in events


@pytest.mark.e2e
def test_get_assay_run_skip_existing(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("assay")
    output_csv = tmp_path / "out" / "assays.csv"
    _ensure_parent(output_csv)
    output_csv.write_text("existing", encoding="utf-8")
    logger_stub = _patch_logger(monkeypatch, get_assay_data)

    def _unexpected_call(*_: object, **__: object) -> int:
        raise AssertionError("run_chembl should be skipped")

    monkeypatch.setattr(get_assay_data, "run_chembl", _unexpected_call)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=True,
        force=False,
    )

    rc = get_assay_data.run(cfg, args)

    assert rc == 0
    events = [event for _, event, _ in logger_stub.events]
    assert "pipeline_skip_existing" in events


@pytest.mark.e2e
def test_get_assay_run_failure(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("assay")
    output_csv = tmp_path / "out" / "assays.csv"
    logger_stub = _patch_logger(monkeypatch, get_assay_data)

    def _failing_run(config: Config, args: argparse.Namespace) -> int:
        get_assay_data.logger.error(
            "assay_pipeline_failed",
            output=str(args.final_out),
            processed=0,
            exit_code=1,
        )
        return 1

    monkeypatch.setattr(get_assay_data, "run_chembl", _failing_run)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
    )

    rc = get_assay_data.run(cfg, args)

    assert rc == 1
    events = [event for _, event, _ in logger_stub.events]
    assert "assay_pipeline_failed" in events


@pytest.mark.e2e
def test_get_activity_run_success(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("activity")
    output_csv = tmp_path / "out" / "activities.csv"
    logger_stub = _patch_logger(monkeypatch, get_activity_data)

    def _stub_run_chembl(config: Config, args: argparse.Namespace) -> int:
        frame = pd.read_csv(args.input_csv)
        frame = frame.drop_duplicates(
            subset=["activity_id", "assay_chembl_id", "molecule_chembl_id"]
        )
        frame["standard_value"] = frame["standard_value"].astype("string").str.strip()
        numeric = pd.to_numeric(frame["standard_value"], errors="coerce")
        missing = numeric.isna()
        if int(missing.sum()):
            get_activity_data.logger.error(
                "activity_missing_value", count=int(missing.sum())
            )
        frame["standard_value_numeric"] = numeric.astype("Float64")
        frame["is_valid"] = (~missing).astype("boolean")
        frame["standard_units"] = frame["standard_units"].astype("string").str.strip()
        frame = frame.sort_values("activity_id").reset_index(drop=True)
        output_path = Path(args.final_out)
        _ensure_parent(output_path)
        frame.to_csv(output_path, index=False)
        get_activity_data.logger.info(
            "activity_pipeline_done", output=str(args.final_out), rows=len(frame)
        )
        return 0

    monkeypatch.setattr(get_activity_data, "run_chembl", _stub_run_chembl)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
    )

    rc = get_activity_data.run(cfg, args)

    assert rc == 0
    result = pd.read_csv(output_csv)
    assert "standard_value_numeric" in result.columns
    assert "is_valid" in result.columns
    assert not result[result["activity_id"] == "A3"]["is_valid"].iloc[0]
    events = [event for _, event, _ in logger_stub.events]
    assert "activity_pipeline_done" in events
    assert "activity_missing_value" in events


@pytest.mark.e2e
def test_get_activity_run_skip_existing(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("activity")
    output_csv = tmp_path / "out" / "activities.csv"
    _ensure_parent(output_csv)
    output_csv.write_text("existing", encoding="utf-8")
    logger_stub = _patch_logger(monkeypatch, get_activity_data)

    def _unexpected_call(*_: object, **__: object) -> int:
        raise AssertionError("run_chembl should be skipped")

    monkeypatch.setattr(get_activity_data, "run_chembl", _unexpected_call)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=True,
        force=False,
    )

    rc = get_activity_data.run(cfg, args)

    assert rc == 0
    events = [event for _, event, _ in logger_stub.events]
    assert "pipeline_skip_existing" in events


@pytest.mark.e2e
def test_get_activity_run_failure(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("activity")
    output_csv = tmp_path / "out" / "activities.csv"
    logger_stub = _patch_logger(monkeypatch, get_activity_data)

    def _failing_run(config: Config, args: argparse.Namespace) -> int:
        get_activity_data.logger.error(
            "activity_pipeline_failed", output=str(args.final_out), exit_code=1
        )
        return 1

    monkeypatch.setattr(get_activity_data, "run_chembl", _failing_run)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
    )

    rc = get_activity_data.run(cfg, args)

    assert rc == 1
    events = [event for _, event, _ in logger_stub.events]
    assert "activity_pipeline_failed" in events


@pytest.mark.e2e
def test_get_activity_run_retry_and_idempotent(
    tmp_path: Path,
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "activities.csv"
    input_csv.write_text("activity_id\nACT1\n", encoding="utf-8")
    output_csv = tmp_path / "out" / "activities.csv"
    logger_stub = _patch_logger(monkeypatch, get_activity_data)

    cfg.activity.dry_run = False
    cfg.activity.limit = None
    cfg.activity.batch_size = 3
    cfg.activity.workers = 1
    cfg.system.doc_quality.enable = False
    cfg.activity_enrichment.action_type.enabled = False
    cfg.activity_enrichment.action_type.log_missing = False
    cfg.activity_enrichment.activity_properties.enabled = False

    frame = pd.DataFrame(
        [
            {
                "activity_id": "ACT1",
                "molecule_chembl_id": "CHEMBL1",
                "assay_chembl_id": "ASSAY1",
                "standard_value": 5.0,
            }
        ]
    )

    call_counter = {"count": 0}

    def _fetch_with_retry(chunk_ids: list[str], **_: object) -> pd.DataFrame:
        call_counter["count"] += 1
        if call_counter["count"] == 1:
            raise requests.RequestException("temporary failure")
        return frame.copy()

    monkeypatch.setattr(get_activity_data.cl, "get_activities", _fetch_with_retry)
    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        _DummyChemblClient,
    )

    written: list[Path] = []

    def _writer_stub(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        *,
        key_cols: Iterable[str] | None,
        col_order: Iterable[str] | None,
        chunksize: int,
        sort_chunksize: int | None = None,
        sep: str = ",",
        encoding: str = "utf-8",
        cfg=None,
        **_: object,
    ) -> Path:
        frames = list(chunks)
        result = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(destination, index=False, sep=sep, encoding=encoding)
        written.append(destination)
        return destination

    monkeypatch.setattr(
        get_activity_data, "write_csv_chunks_deterministic", _writer_stub
    )

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
    )

    rc = get_activity_data.run(cfg, args)

    assert rc == 0
    assert call_counter["count"] == 2
    events = [event for _, event, _ in logger_stub.events]
    assert "activity_fetch_retry" in events
    assert "activity_pipeline_done" in events
    assert output_csv in written

    args.skip_existing = True
    rc_second = get_activity_data.run(cfg, args)

    assert rc_second == 0
    events_second = [event for _, event, _ in logger_stub.events]
    assert "pipeline_skip_existing" in events_second


@pytest.mark.e2e
def test_get_activity_run_workers_offset_and_non_csv(
    tmp_path: Path,
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "activities.csv"
    input_csv.write_text("activity_id\nACT1\nACT2\n", encoding="utf-8")
    output_csv = tmp_path / "out" / "activities.tsv"
    logger_stub = _patch_logger(monkeypatch, get_activity_data)

    cfg.activity.dry_run = False
    cfg.activity.limit = None
    cfg.activity.batch_size = 2
    cfg.activity.workers = 3
    cfg.system.doc_quality.enable = False
    cfg.activity_enrichment.action_type.enabled = False
    cfg.activity_enrichment.action_type.log_missing = False
    cfg.activity_enrichment.activity_properties.enabled = False
    cfg.io.csv_sep = "\t"

    monkeypatch.setattr(
        get_activity_data.io,
        "read_ids",
        lambda *_args, **_kwargs: iter(["ACT0", "ACT1", "ACT2"]),
    )

    captured: dict[str, int] = {}

    def _fake_pipeline(**kwargs):
        fetch_config = kwargs["fetch_config"]
        captured["workers"] = fetch_config.workers
        output_path = kwargs["output_path"]
        output_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            [
                {
                    "activity_id": "ACT2",
                    "molecule_chembl_id": "CHEMBL2",
                    "assay_chembl_id": "ASSAY2",
                    "standard_value": 7.0,
                }
            ]
        ).to_csv(output_path, index=False, sep=cfg.io.csv_sep)
        return PipelineRunResult(exit_code=0, output_path=output_path, written=True)

    monkeypatch.setattr(
        "library.pipelines.activity.run.run_activity_pipeline",
        _fake_pipeline,
    )
    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        _DummyChemblClient,
    )

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
        offset=1,
        workers=cfg.activity.workers,
    )

    rc = get_activity_data.run(cfg, args)

    assert rc == 0
    events = [event for _, event, _ in logger_stub.events]
    assert "activity_pipeline_done" in events
    assert any(
        event_name == "process_offset" for _, event_name, _ in logger_stub.events
    )
    assert captured["workers"] == max(1, cfg.activity.workers)
    assert output_csv.exists()
