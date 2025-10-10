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
def test_resolve_chunk_size__cases(
    value: int | None, expected: int | None, logger_stub: _MemoryLogger
) -> None:
    result = get_document_data._resolve_chunk_size(value)

    assert result == expected
    if value is not None and value <= 0:
        assert (
            "warning",
            "invalid_csv_chunksize",
            {"value": value},
        ) in logger_stub.events


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


def test_run__skip_existing(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
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
    assert (
        "info",
        "pipeline_skip_existing",
        {"output": str(output_csv)},
    ) in logger_stub.events


def test_run__missing_handler_logs_error(
    cfg: Config, tmp_path: Path, logger_stub: _MemoryLogger
) -> None:
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


def test_run__propagates_timeout(
    cfg: Config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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


def test_run__pubmed_timeout_override_updates_config(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("PMID\n123\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    recorded: dict[str, float] = {}

    def fake_pubmed(
        cfg_obj: Config, args: argparse.Namespace, *, pipeline=None
    ) -> int:  # noqa: ARG001
        recorded["pubmed"] = cfg_obj.pubmed.timeout_read
        recorded["api"] = cfg_obj.api.timeout_read
        return 0

    monkeypatch.setitem(get_document_data.MODE_HANDLERS, "pubmed", fake_pubmed)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
        command="pubmed",
        mode="pubmed",
        func=fake_pubmed,
        timeout=25.0,
        fallback_doi_enabled=False,
        fallback_doi_overwrite=False,
        fallback_doi_path=None,
    )

    previous_api_timeout = cfg.api.timeout_read

    exit_code = get_document_data.run(cfg, args)

    assert exit_code == 0
    assert recorded["pubmed"] == pytest.approx(25.0)
    assert recorded["api"] == pytest.approx(previous_api_timeout)


def test_run__all_mode_pubmed_timeout_option_updates_config(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\n1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    recorded: dict[str, float] = {}

    def fake_all(
        cfg_obj: Config, args: argparse.Namespace, *, pipeline=None
    ) -> int:  # noqa: ARG001
        recorded["pubmed"] = cfg_obj.pubmed.timeout_read
        recorded["api"] = cfg_obj.api.timeout_read
        recorded["pubmed_arg"] = args.pubmed_timeout
        recorded["chembl_arg"] = getattr(args, "chembl_timeout", None)
        return 0

    monkeypatch.setitem(get_document_data.MODE_HANDLERS, "all", fake_all)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
        command="all",
        mode="all",
        func=fake_all,
        timeout=None,
        pubmed_timeout=33.0,
        chembl_timeout=None,
        fallback_doi_enabled=False,
        fallback_doi_overwrite=False,
        fallback_doi_path=None,
    )

    exit_code = get_document_data.run(cfg, args)

    assert exit_code == 0
    assert recorded["pubmed"] == pytest.approx(33.0)
    assert recorded["pubmed_arg"] == pytest.approx(33.0)
    assert recorded["api"] == pytest.approx(cfg.api.timeout_read)


def test_run__all_mode_timeout_fallback_updates_pubmed(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\n1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    recorded: dict[str, float] = {}

    def fake_all(
        cfg_obj: Config, args: argparse.Namespace, *, pipeline=None
    ) -> int:  # noqa: ARG001
        recorded["pubmed"] = cfg_obj.pubmed.timeout_read
        recorded["pubmed_arg"] = args.pubmed_timeout
        return 0

    monkeypatch.setitem(get_document_data.MODE_HANDLERS, "all", fake_all)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
        command="all",
        mode="all",
        func=fake_all,
        timeout=41.0,
        pubmed_timeout=None,
        chembl_timeout=None,
        fallback_doi_enabled=False,
        fallback_doi_overwrite=False,
        fallback_doi_path=None,
    )

    exit_code = get_document_data.run(cfg, args)

    assert exit_code == 0
    assert recorded["pubmed"] == pytest.approx(41.0)
    value = recorded["pubmed_arg"]
    if value is not None:
        assert value == pytest.approx(41.0)


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

    def fail_postprocessing(
        path: Path, *, skip_qa: bool = False
    ) -> None:  # noqa: ARG001
        assert skip_qa is False
        raise RuntimeError("QA mismatches")

    monkeypatch.setattr(
        get_document_data,
        "_maybe_run_document_postprocessing",
        fail_postprocessing,
    )

    result = get_document_data._finalise_export(
        frame,
        output_csv,
        cfg,
        input_csv=input_csv,
    )

    assert result.exit_code == 1
    assert (
        "error",
        "document_postprocess_qa_mismatch",
        {"error": "QA mismatches", "path": str(output_csv)},
    ) in logger_stub.events
    assert finalise_calls, "finalise_csv_output should still be invoked"


def test_finalise_export__partial_run_skips_qa(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    data_dir = tmp_path / "data"
    reference_dir = data_dir / "input" / "full"
    reference_dir.mkdir(parents=True, exist_ok=True)
    (reference_dir / "document.csv").write_text("reference", encoding="utf-8")

    output_csv = data_dir / "output" / "document" / "output.document_20250101.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("document_chembl_id\nCHEMBL1\n", encoding="utf-8")

    frame = get_document_data.build_dataframe(
        [{"document_chembl_id": "CHEMBL1"}],
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

    monkeypatch.setattr(
        get_document_data,
        "finalise_csv_output",
        lambda **_: None,
    )

    captured: dict[str, Any] = {}

    def fake_preprocess(
        base_path: str,
        ref_document_rel: str = "",
        out_document_rel: str = "",
        qa_reference_rel: str | None = None,
        *,
        run_qa: bool = True,
    ) -> str:  # noqa: ARG001
        captured.update(
            {
                "base_path": base_path,
                "ref_rel": ref_document_rel,
                "out_rel": out_document_rel,
                "run_qa": run_qa,
            }
        )
        return str(Path(base_path) / "output" / "document" / "preprocessed.csv")

    monkeypatch.setattr(
        get_document_data,
        "preprocess_documents_csv",
        fake_preprocess,
    )

    result = get_document_data._finalise_export(
        frame,
        output_csv,
        cfg,
        input_csv=input_csv,
        key_columns=["document_chembl_id"],
        partial_run=True,
    )

    assert result.exit_code == 0
    assert captured.get("run_qa") is False
    assert captured.get("base_path") == str(data_dir)
    assert (
        "info",
        "document_postprocess_qa_skipped_partial",
        {"output": str(output_csv), "reason": "partial_run"},
    ) in logger_stub.events


def test_finalise_export__reuses_postprocess_when_rerun_disabled(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_csv = tmp_path / "output.document_20250101.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("document_chembl_id\nCHEMBL1\n", encoding="utf-8")

    frame = get_document_data.build_dataframe(
        [{"document_chembl_id": "CHEMBL1"}],
        columns=get_document_data.DOCUMENT_SCHEMA_COLUMNS,
        fill_missing=False,
    )

    def fake_write_csv_chunks(
        chunks: Iterable[pd.DataFrame],
        path: Path,
        **kwargs: Any,
    ) -> Path:
        frames = list(chunks)
        combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        resolved = Path(path)
        resolved.parent.mkdir(parents=True, exist_ok=True)
        combined.to_csv(resolved, index=False)
        return resolved

    monkeypatch.setattr(
        get_document_data,
        "write_csv_chunks_deterministic",
        fake_write_csv_chunks,
    )
    monkeypatch.setattr(get_document_data, "finalise_csv_output", lambda **_: None)

    expected_post = get_document_data.document_export_postprocessing.resolve_default_destination(
        output_csv
    )
    expected_post.parent.mkdir(parents=True, exist_ok=True)
    expected_post.write_text("header\n", encoding="utf-8")

    call_count = 0

    def _unexpected_postprocess(path: Path, *, cfg: Any) -> Path:  # noqa: ARG001
        nonlocal call_count
        call_count += 1
        return Path(path)

    monkeypatch.setattr(
        get_document_data.document_export_postprocessing,
        "postprocess_export_file",
        _unexpected_postprocess,
    )

    result = get_document_data._finalise_export(
        frame,
        output_csv,
        cfg,
        input_csv=input_csv,
        rerun_postprocess=False,
    )

    assert result.exit_code == 0
    assert result.postprocess_path == expected_post
    assert call_count == 0
    assert (
        "info",
        "document_export_postprocess_written",
        {"path": str(expected_post)},
    ) in logger_stub.events


def test_finalise_export__rerun_enabled_invokes_postprocess(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_csv = tmp_path / "output.document_20250101.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("document_chembl_id\nCHEMBL1\n", encoding="utf-8")

    frame = get_document_data.build_dataframe(
        [{"document_chembl_id": "CHEMBL1"}],
        columns=get_document_data.DOCUMENT_SCHEMA_COLUMNS,
        fill_missing=False,
    )

    def fake_write_csv_chunks(
        chunks: Iterable[pd.DataFrame],
        path: Path,
        **kwargs: Any,
    ) -> Path:
        frames = list(chunks)
        combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        resolved = Path(path)
        resolved.parent.mkdir(parents=True, exist_ok=True)
        combined.to_csv(resolved, index=False)
        return resolved

    monkeypatch.setattr(
        get_document_data,
        "write_csv_chunks_deterministic",
        fake_write_csv_chunks,
    )
    monkeypatch.setattr(get_document_data, "finalise_csv_output", lambda **_: None)

    generated_post = output_csv.with_name("preprocessed_output.document_20250101.csv")

    def fake_postprocess(path: Path, *, cfg: Any) -> Path:  # noqa: ARG001
        generated_post.write_text("header\n", encoding="utf-8")
        return generated_post

    call_count = 0

    def counting_postprocess(path: Path, *, cfg: Any) -> Path:  # noqa: ARG001
        nonlocal call_count
        call_count += 1
        return fake_postprocess(path, cfg=cfg)

    monkeypatch.setattr(
        get_document_data.document_export_postprocessing,
        "postprocess_export_file",
        counting_postprocess,
    )

    # Create an existing file to confirm rerun still executes the postprocess stage.
    generated_post.parent.mkdir(parents=True, exist_ok=True)
    generated_post.write_text("old\n", encoding="utf-8")

    result = get_document_data._finalise_export(
        frame,
        output_csv,
        cfg,
        input_csv=input_csv,
        rerun_postprocess=True,
    )

    assert result.exit_code == 0
    assert result.postprocess_path == generated_post
    assert call_count == 1
    assert (
        "info",
        "document_export_postprocess_written",
        {"path": str(generated_post)},
    ) in logger_stub.events
