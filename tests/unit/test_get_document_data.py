"""Unit tests for :mod:`scripts.get_document_data`."""

from __future__ import annotations

import argparse
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library.config import Config
from scripts import get_document_data


class _MemoryLogger:
    """Capture structured log events for assertions."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def info(self, event: str, **payload: object) -> None:
        self.events.append(("info", event, dict(payload)))

    def warning(self, event: str, **payload: object) -> None:
        self.events.append(("warning", event, dict(payload)))

    def error(self, event: str, **payload: object) -> None:
        self.events.append(("error", event, dict(payload)))


@pytest.fixture()
def logger_stub(monkeypatch: pytest.MonkeyPatch) -> _MemoryLogger:
    logger = _MemoryLogger()
    monkeypatch.setattr(get_document_data, "logger", logger)
    return logger


@pytest.mark.parametrize(
    "value, expected, event",
    [
        (None, None, None),
        (5, 5, None),
        (True, None, "invalid_stream_chunk_size_bool"),
        (" 7 ", 7, None),
        ("", None, None),
        ("abc", None, "invalid_stream_chunk_size_string"),
        (2.0, 2, None),
        (2.5, None, "invalid_stream_chunk_size_float"),
        (float("nan"), None, None),
        (pd.Series([3]), 3, None),
        (pd.Series([1, 2]), None, "invalid_stream_chunk_size_series"),
        ({"size": 1}, None, "invalid_stream_chunk_size_type"),
    ],
)
def test_coerce_chunk_size_value__cases(
    value: object,
    expected: int | None,
    event: str | None,
    logger_stub: _MemoryLogger,
) -> None:
    result = get_document_data._coerce_chunk_size_value(value)

    assert result == expected
    if event is not None:
        assert any(evt == event for _, evt, _ in logger_stub.events)


@pytest.mark.parametrize(
    "value, expected",
    [
        (None, None),
        (0, None),
        (10, 10),
        (-5, None),
    ],
)
def test_resolve_chunk_size__cases(value: int | None, expected: int | None, logger_stub: _MemoryLogger) -> None:
    result = get_document_data._resolve_chunk_size(value)

    assert result == expected
    if value is not None and value <= 0:
        assert ("warning", "invalid_csv_chunksize", {"value": value}) in logger_stub.events


def test_coalesce_columns__prefers_first_non_empty() -> None:
    frame = pd.DataFrame(
        {
            "first": ["", "primary", None],
            "second": ["fallback", "", "replacement"],
            "third": ["extra", "extra", ""],
        }
    )

    result = get_document_data._coalesce_columns(frame, ["first", "second", "third"])

    assert result.tolist() == ["fallback", "primary", "replacement"]


def test_resolve_duplicate_column__merges_frames() -> None:
    frame = pd.DataFrame(
        [["first", "", "x"], ["", "second", "y"]],
        columns=["title", "title", "other"],
    )

    result = get_document_data._resolve_duplicate_column(frame, "title")

    assert result.tolist() == ["first", "second"]


def test_collapse_duplicate_columns__projects_first_instance() -> None:
    frame = pd.DataFrame(
        [
            ["1", "1", "10.1", "Primary", "Primary"],
            ["2", "2", "", "Secondary", "Secondary"],
        ],
        columns=[
            "PubMed.PMID",
            "PubMed.PMID",
            "PubMed.DOI",
            "title",
            "title",
        ],
    )

    result = get_document_data._collapse_duplicate_columns(frame)

    assert list(result.columns) == ["PubMed.PMID", "PubMed.DOI", "title"]
    assert result["title"].tolist() == ["Primary", "Secondary"]


def test_prepare_export_frame__renames_and_coalesces() -> None:
    frame = pd.DataFrame(
        {
            "ChEMBL.document_chembl_id": ["CHEMBL1"],
            "ChEMBL.pubmed_id": ["123"],
            "ChEMBL.doi": ["10.1/abc"],
            "title": [" Example "],
            "abstract": ["Test"],
        }
    )

    export = get_document_data._prepare_export_frame(frame)

    assert "PubMed.PMID" in export.columns
    assert export.loc[0, "PubMed.PMID"] == ""
    assert export.loc[0, "ChEMBL.pubmed_id"] == "123"
    assert export.loc[0, "ChEMBL.doi"] == "10.1/abc"
    assert export.loc[0, "ChEMBL.title"] == " Example "


def test_run__skip_existing(cfg: Config, tmp_path: Path, logger_stub: _MemoryLogger, monkeypatch: pytest.MonkeyPatch) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\n1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"
    output_csv.write_text("existing", encoding="utf-8")

    called = False

    def fake_run(cfg: Config, args: argparse.Namespace) -> int:
        nonlocal called
        called = True
        return 0

    monkeypatch.setattr(get_document_data, "run_all", fake_run)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=True,
        force=False,
        command="all",
        func=fake_run,
        timeout=None,
    )

    exit_code = get_document_data.run(cfg, args)

    assert exit_code == 0
    assert not called
    assert ("info", "pipeline_skip_existing", {"output": str(output_csv)}) in logger_stub.events


def test_run__missing_handler_logs_error(cfg: Config, tmp_path: Path, logger_stub: _MemoryLogger) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\n1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
        command="all",
        func=None,
        timeout=None,
    )

    exit_code = get_document_data.run(cfg, args)

    assert exit_code == 1
    assert (
        "error",
        "missing_subcommand_handler",
        {"command": "all"},
    ) in logger_stub.events


def test_run__propagates_timeout(cfg: Config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\n1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    called: list[float | None] = []

    def fake_run(cfg: Config, args: argparse.Namespace) -> int:
        called.append(cfg.api.timeout_read)
        return 0

    monkeypatch.setattr(get_document_data, "run_all", fake_run)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
        command="all",
        func=fake_run,
        timeout=42.5,
    )

    exit_code = get_document_data.run(cfg, args)

    assert exit_code == 0
    assert called == [42.5]


def test_finalise_export__qa_mismatch_sets_exit_code(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_csv = tmp_path / "output.document_20250101.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("document_chembl_id\n", encoding="utf-8")

    frame = get_document_data.build_dataframe(
        [],
        columns=get_document_data.DOCUMENT_SCHEMA_COLUMNS,
        fill_missing=False,
    )

    def fake_write_csv_chunks(
        chunks: Iterable[pd.DataFrame],
        path: Path,
        **kwargs: Any,
    ) -> Path:
        frames = list(chunks)
        if frames:
            combined = pd.concat(frames, ignore_index=True)
        else:
            combined = pd.DataFrame(columns=kwargs.get("col_order", []))
        resolved = Path(path)
        resolved.parent.mkdir(parents=True, exist_ok=True)
        combined.to_csv(resolved, index=False)
        return resolved

    monkeypatch.setattr(
        get_document_data,
        "write_csv_chunks_deterministic",
        fake_write_csv_chunks,
    )

    def fake_postprocess_export(path: Path, *, cfg: Any) -> Path:  # noqa: ARG001
        return Path(path)

    monkeypatch.setattr(
        get_document_data.document_export_postprocessing,
        "postprocess_export_file",
        fake_postprocess_export,
    )

    finalise_calls: list[dict[str, Any]] = []

    def fake_finalise_csv_output(**kwargs: Any) -> None:
        finalise_calls.append(kwargs)

    monkeypatch.setattr(
        get_document_data,
        "finalise_csv_output",
        fake_finalise_csv_output,
    )

    def fail_postprocessing(path: Path) -> None:  # noqa: ARG001
        raise RuntimeError("QA mismatches")

    monkeypatch.setattr(
        get_document_data,
        "_maybe_run_document_postprocessing",
        fail_postprocessing,
    )

    exit_code = get_document_data._finalise_export(
        frame,
        output_csv,
        cfg,
        input_csv=input_csv,
    )

    assert exit_code == 1
    assert (
        "error",
        "document_postprocess_qa_mismatch",
        {"error": "QA mismatches", "path": str(output_csv)},
    ) in logger_stub.events
    assert finalise_calls, "finalise_csv_output should still be invoked"
