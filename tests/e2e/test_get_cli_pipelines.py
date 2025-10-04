"""Integration tests for the CLI entrypoints in ``scripts/get_*`` modules."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from typing import Callable

import pandas as pd
import pytest

from library.config import Config
from scripts import (
    get_activity_data,
    get_assay_data,
    get_document_data,
    get_target_data,
    get_testitem_data,
)


@pytest.fixture()
def sample_csv(tmp_path: Path) -> Callable[[str], Path]:
    """Return a helper copying CSV fixtures into ``tmp_path``."""

    data_dir = Path(__file__).resolve().parents[1] / "data"

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

    def info(self, event: str, **kwargs: object) -> None:
        self.events.append(("info", event, dict(kwargs)))

    def warning(self, event: str, **kwargs: object) -> None:
        self.events.append(("warning", event, dict(kwargs)))

    def error(self, event: str, **kwargs: object) -> None:
        self.events.append(("error", event, dict(kwargs)))


def _patch_logger(monkeypatch: pytest.MonkeyPatch, module: object) -> _MemoryLogger:
    logger = _MemoryLogger()
    monkeypatch.setattr(module, "logger", logger)
    return logger

@pytest.mark.e2e
def test_get_testitem_run_success(
    tmp_path: Path,
    sample_csv: Callable[[str], Path],
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = sample_csv("testitem")
    output_csv = tmp_path / "out" / "testitems.csv"
    logger_stub = _patch_logger(monkeypatch, get_testitem_data)

    def _stub_pipeline(config: Config, options: get_testitem_data.TestitemPipelineOptions) -> int:
        assert config is cfg
        assert options.input_csv == Path(input_csv)
        assert options.output_csv == output_csv
        frame = pd.read_csv(options.input_csv)
        frame["preferred_name"] = frame["preferred_name"].fillna("").astype("string").str.strip()
        missing = frame["preferred_name"] == ""
        if int(missing.sum()):
            get_testitem_data.logger.warning(
                "testitem_missing_name", count=int(missing.sum())
            )
        frame["normalized_name"] = frame["preferred_name"].replace("", pd.NA).str.lower()
        frame["is_named"] = (~missing).astype("boolean")
        _ensure_parent(output_csv)
        frame.to_csv(output_csv, index=False)
        return 0

    monkeypatch.setattr(get_testitem_data, "run_testitem_pipeline", _stub_pipeline)

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        skip_existing=False,
        force=False,
    )

    rc = get_testitem_data.run(cfg, args)

    assert rc == 0
    result = pd.read_csv(output_csv)
    assert list(result.columns) == [
        "molecule_chembl_id",
        "preferred_name",
        "normalized_name",
        "is_named",
    ]
    assert result.loc[0, "normalized_name"] == "compound 1"
    assert pd.isna(result.loc[1, "normalized_name"])
    assert not result.loc[1, "is_named"]
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
    output_csv = tmp_path / "out" / "testitems.csv"
    logger_stub = _patch_logger(monkeypatch, get_testitem_data)

    def _failing_pipeline(config: Config, options: get_testitem_data.TestitemPipelineOptions) -> int:
        get_testitem_data.logger.error(
            "testitem_pipeline_failed", output=str(options.output_csv), exit_code=2
        )
        return 2

    monkeypatch.setattr(get_testitem_data, "run_testitem_pipeline", _failing_pipeline)

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
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
    output_csv = tmp_path / "out" / "testitems.csv"
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
        output_csv=output_csv,
        skip_existing=True,
        force=False,
    )

    rc = get_testitem_data.run(cfg, args)

    assert rc == 0
    assert call_counter["called"] == 0
    events = [event for _, event, _ in logger_stub.events]
    assert "pipeline_skip_existing" in events


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
        _ensure_parent(Path(args.output_csv))
        frame.to_csv(args.output_csv, index=False)
        get_document_data.logger.info(
            "document_all_done", output=str(args.output_csv), rows=len(frame)
        )
        return 0

    monkeypatch.setattr(get_document_data, "run_all", _stub_all)

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
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
    assert not result.loc[result["document_chembl_id"] == "CHEMBL456", "has_pubmed"].iloc[0]
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
        output_csv=output_csv,
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
            "document_all_failed", output=str(args.output_csv), exit_code=1
        )
        return 1

    monkeypatch.setattr(get_document_data, "run_all", _failing_all)

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
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

    def _stub_run_chembl(config: Config, args: argparse.Namespace) -> int:
        frame = pd.read_csv(args.input_csv)
        frame["description"] = frame["description"].astype("string").str.strip()
        frame["description_length"] = frame["description"].str.len().astype("Int64")
        frame = frame.drop_duplicates(subset=["assay_chembl_id"])
        frame = frame.sort_values("assay_chembl_id").reset_index(drop=True)
        _ensure_parent(args.output_csv)
        frame.to_csv(args.output_csv, index=False)
        get_assay_data.logger.info(
            "assay_pipeline_done", output=str(args.output_csv), processed=len(frame)
        )
        return 0

    monkeypatch.setattr(get_assay_data, "run_chembl", _stub_run_chembl)

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
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
    ]
    assert (result["description_length"] == result["description"].str.len()).all()
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
        output_csv=output_csv,
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
            "assay_pipeline_failed", output=str(args.output_csv), processed=0, exit_code=1
        )
        return 1

    monkeypatch.setattr(get_assay_data, "run_chembl", _failing_run)

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
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
        frame = frame.drop_duplicates(subset=["activity_id", "assay_chembl_id", "molecule_chembl_id"])
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
        _ensure_parent(args.output_csv)
        frame.to_csv(args.output_csv, index=False)
        get_activity_data.logger.info(
            "activity_pipeline_done", output=str(args.output_csv), rows=len(frame)
        )
        return 0

    monkeypatch.setattr(get_activity_data, "run_chembl", _stub_run_chembl)

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
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
        output_csv=output_csv,
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
            "activity_pipeline_failed", output=str(args.output_csv), exit_code=1
        )
        return 1

    monkeypatch.setattr(get_activity_data, "run_chembl", _failing_run)

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        skip_existing=False,
        force=False,
    )

    rc = get_activity_data.run(cfg, args)

    assert rc == 1
    events = [event for _, event, _ in logger_stub.events]
    assert "activity_pipeline_failed" in events
