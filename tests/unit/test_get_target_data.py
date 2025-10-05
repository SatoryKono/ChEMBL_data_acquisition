"""Unit tests for :mod:`scripts.get_target_data`."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd
from types import SimpleNamespace

import pytest

from library.config import Config
from scripts import get_target_data


class _MemoryLogger:
    """Capture structured log events emitted by the target pipeline."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def info(self, event: str, **payload: object) -> None:
        self.events.append(("info", event, dict(payload)))

    def warning(self, event: str, **payload: object) -> None:
        self.events.append(("warning", event, dict(payload)))

    def error(self, event: str, **payload: object) -> None:
        self.events.append(("error", event, dict(payload)))

    def exception(self, event: str, **payload: object) -> None:
        self.events.append(("exception", event, dict(payload)))


@pytest.fixture()
def logger_stub(monkeypatch: pytest.MonkeyPatch) -> _MemoryLogger:
    logger = _MemoryLogger()
    monkeypatch.setattr(get_target_data, "logger", logger)
    return logger


@pytest.mark.parametrize(
    "latin, cyrillic",
    list(get_target_data._RUSSIAN_KEYBOARD_MAP.items())[:10],
)
def test_translate_keyboard_layout__single_characters(latin: str, cyrillic: str) -> None:
    assert get_target_data._translate_keyboard_layout(latin) == cyrillic
    assert get_target_data._translate_keyboard_layout(latin.upper()) == cyrillic.upper()


@pytest.mark.parametrize("text", ["Hello", "Test-Value", "123", "aA" ])
def test_translate_keyboard_layout__composite_strings(text: str) -> None:
    translated = get_target_data._translate_keyboard_layout(text)
    assert len(translated) == len(text)
    for original, result in zip(text, translated):
        lower = original.lower()
        mapped = get_target_data._RUSSIAN_KEYBOARD_MAP.get(lower, lower)
        expected = mapped.upper() if original.isupper() else mapped
        assert result == expected


@pytest.mark.parametrize("command", ["fetch", "all", "run"])
def test_keyboard_aliases__cases(command: str) -> None:
    aliases = get_target_data._keyboard_aliases(command)
    translated = get_target_data._translate_keyboard_layout(command)
    expected_variants = {translated, translated.capitalize(), translated.upper()}
    expected_variants.discard(command)
    assert set(aliases) == expected_variants


@pytest.mark.parametrize(
    "value, tokens",
    [
        ("P12345", ["P12345"]),
        ("P12345|Q67890", ["P12345", "Q67890"]),
        (" |P99999| ", ["P99999"]),
        ("-", []),
        ("P12345-2|-|Q99999", ["P12345-2", "Q99999"]),
    ],
)
def test_split_uniprot_tokens__cases(value: str, tokens: list[str]) -> None:
    assert list(get_target_data._split_uniprot_tokens(value)) == tokens


def test_collect_uniprot_candidate_columns__orders_columns(cfg: Config) -> None:
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            cfg.target.all.uniprot_column: ["P12345"],
            "mapping_uniprot_id": ["Q99999"],
            "extra_accession": ["E11111"],
            "uniprot_secondary": ["S22222"],
        }
    )

    ordered = get_target_data._collect_uniprot_candidate_columns(frame, cfg)

    assert ordered[0] == cfg.target.all.uniprot_column
    assert "mapping_uniprot_id" in ordered
    assert any("accession" in column for column in ordered)


def test_ensure_merge_column_present__uses_alias(cfg: Config, logger_stub: _MemoryLogger) -> None:
    cfg.target.all.uniprot_column = "canonical_uniprot"
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "mapping_uniprot_id": ["P12345"],
        }
    )

    result = get_target_data._ensure_merge_column_present(frame, cfg.target.all.uniprot_column, cfg)

    assert "canonical_uniprot" in result.columns
    assert result.loc[0, "canonical_uniprot"] == "P12345"
    assert any(event == "uniprot_merge_column_alias" for _, event, _ in logger_stub.events)


def test_ensure_merge_column_present__raises(cfg: Config, logger_stub: _MemoryLogger) -> None:
    cfg.target.all.uniprot_column = "canonical_uniprot"
    frame = pd.DataFrame({"target_chembl_id": ["CHEMBL1"]})

    with pytest.raises(get_target_data.PipelineError):
        get_target_data._ensure_merge_column_present(frame, cfg.target.all.uniprot_column, cfg)

    assert any(event == "missing_uniprot_column" for _, event, _ in logger_stub.events)


@pytest.mark.parametrize(
    "values, expected",
    [
        (("P12345", "Q67890"), "P12345|Q67890"),
        (("", "Q67890"), "Q67890"),
        ((None, ""), ""),
        (("A", None, "B"), "A|B"),
    ],
)
def test_pipe_merge__cases(values: Iterable[str | None], expected: str) -> None:
    assert get_target_data._pipe_merge(list(values)) == expected


@pytest.mark.parametrize(
    "value, expected",
    [
        ("P12345", "P12345"),
        ("P12345|Q67890", "P12345"),
        (None, ""),
        (" ", " "),
    ],
)
def test_first_token__cases(value: str | None, expected: str) -> None:
    assert get_target_data._first_token(value) == expected


def test_limited_ids__respects_limit() -> None:
    source = iter(["A", "B", "C"])

    limited = list(get_target_data._limited_ids(source, 2))

    assert limited == ["A", "B"]


def test_run__skip_existing(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\n1\n", encoding="utf-8")
    final_out = tmp_path / "targets.csv"
    final_out.write_text("existing", encoding="utf-8")

    called = False

    def fake_run(cfg: Config, args: argparse.Namespace) -> int:
        nonlocal called
        called = True
        return 0

    monkeypatch.setattr(get_target_data, "run_all", fake_run)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=final_out,
        raw_out=None,
        raw_format="csv",
        no_reindex_raw=False,
        skip_existing=True,
        force=False,
        command="all",
        func=fake_run,
    )

    exit_code = get_target_data.run(cfg, args)

    assert exit_code == 0
    assert not called
    assert ("info", "pipeline_skip_existing", {"output": str(final_out)}) in logger_stub.events


def test_run_uniprot__invokes_target_postprocess(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg.system.doc_quality.enable = False
    input_csv = tmp_path / "uniprot_ids.csv"
    input_csv.write_text("uniprot_id\nP12345\n", encoding="utf-8")
    output_csv = tmp_path / "output.target_20250101.csv"

    monkeypatch.setattr(get_target_data.uu, "init_session", lambda *_, **__: None)

    def _fake_process(**kwargs: object) -> None:
        assert kwargs["output_csv"] == str(output_csv)
        Path(kwargs["output_csv"]).write_text(
            "uniprot_id\nP12345\n", encoding=cfg.io.csv_encoding
        )

    monkeypatch.setattr(get_target_data.uu, "process", _fake_process)

    recorded: dict[str, object] = {}

    def _fake_postprocess(
        path: Path,
        *,
        cfg: Config,
        context: get_target_data.IsoformPostprocessContext | None = None,
        ambiguous_classifications: int | None = None,
    ) -> None:
        recorded["path"] = Path(path)
        recorded["context"] = context
        recorded["ambiguous"] = ambiguous_classifications

    monkeypatch.setattr(
        get_target_data,
        "_postprocess_target_exports",
        _fake_postprocess,
    )

    args = argparse.Namespace(input_csv=input_csv, output_csv=output_csv)

    exit_code = get_target_data.run_uniprot(cfg, args)

    assert exit_code == 0
    assert recorded["path"] == output_csv
    assert isinstance(recorded["context"], get_target_data.IsoformPostprocessContext)
    assert recorded["context"].args is args
    assert recorded["ambiguous"] is None


def test_run__delegates_to_handler(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\n1\n", encoding="utf-8")
    final_out = tmp_path / "targets.csv"

    monkeypatch.setattr(get_target_data, "run_all", lambda *_: 4)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=final_out,
        raw_out=None,
        raw_format="csv",
        no_reindex_raw=False,
        skip_existing=False,
        force=False,
        command="all",
        func=get_target_data.run_all,
    )

    exit_code = get_target_data.run(cfg, args)

    assert exit_code == 4


def test_postprocess_organism_export__success(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "output.target_20250101.csv"
    source.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    output_path = tmp_path / "organism.output.target_20250101.csv"

    def _fake_postprocess(path: str) -> str:
        assert path == str(source)
        output_path.write_text("target_chembl_id,target_type\nCHEMBL1,Multicellular\n", encoding="utf-8")
        return str(output_path)

    monkeypatch.setattr(get_target_data.target_pp, "postprocess_target_table", _fake_postprocess)

    result = get_target_data._postprocess_organism_export(source, cfg=cfg)

    assert result == output_path
    assert (
        "info",
        "target_organism_postprocess_done",
        {"path": str(output_path), "source": str(source)},
    ) in logger_stub.events


def test_postprocess_organism_export__failure(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "output.target_20250101.csv"
    source.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    def _failing_postprocess(path: str) -> str:
        raise RuntimeError("boom")

    monkeypatch.setattr(get_target_data.target_pp, "postprocess_target_table", _failing_postprocess)

    result = get_target_data._postprocess_organism_export(source, cfg=cfg)

    assert result is None
    assert (
        "exception",
        "target_organism_postprocess_failed",
        {"path": str(source), "error": "boom"},
    ) in logger_stub.events


def test_postprocess_isoform_export__skips_for_custom_name(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "custom_export.csv"
    source.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    def _unexpected(*_: object, **__: object) -> None:  # pragma: no cover - defensive
        raise AssertionError("isoform post-processing should be skipped")

    monkeypatch.setattr(get_target_data.target_pp, "process_targets", _unexpected)

    result = get_target_data._postprocess_isoform_export(source, cfg=cfg)

    assert result is None
    assert (
        "info",
        "target_isoform_postprocess_skipped",
        {"path": str(source), "reason": "unsupported_export_name"},
    ) in logger_stub.events


def test_postprocess_target_exports__chains_helpers(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "output.target_20250101.csv"
    source.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    call_order: list[str] = []

    def _fake_organism(path: Path, *, cfg: Config) -> Path:
        assert Path(path) == source
        call_order.append("organism")
        return source

    def _fake_isoform(
        path: Path,
        *,
        cfg: Config,
        context: get_target_data.IsoformPostprocessContext | None = None,
        ambiguous_classifications: int | None = None,
    ) -> Path | None:
        assert Path(path) == source
        call_order.append("isoform")
        return source

    monkeypatch.setattr(get_target_data, "_postprocess_organism_export", _fake_organism)
    monkeypatch.setattr(get_target_data, "_postprocess_isoform_export", _fake_isoform)

    get_target_data._postprocess_target_exports(source, cfg=cfg)

    assert call_order == ["organism", "isoform"]


def test_postprocess_target_exports__skips_for_uniprot_suffix(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    logger_stub: _MemoryLogger,
) -> None:
    source = tmp_path / "output.targets_20250101_uniprot.csv"
    source.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    def _unexpected(*_: object, **__: object) -> None:  # pragma: no cover - defensive
        raise AssertionError("post-processing should be skipped")

    monkeypatch.setattr(get_target_data, "_postprocess_organism_export", _unexpected)
    monkeypatch.setattr(get_target_data, "_postprocess_isoform_export", _unexpected)
    monkeypatch.setattr(get_target_data, "_postprocess_names_export", _unexpected)
    monkeypatch.setattr(get_target_data, "_postprocess_iuphar_export", _unexpected)

    get_target_data._postprocess_target_exports(source, cfg=cfg)

    assert (
        "info",
        "target_postprocess_skipped",
        {"path": str(source), "reason": "unsupported_export_name"},
    ) in logger_stub.events


def test_postprocess_iuphar_export__missing_columns_logs_warning(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, logger_stub: _MemoryLogger
) -> None:
    source = tmp_path / "output.target_20250101.csv"
    source.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    class _StubError(RuntimeError):
        pass

    def _raise_error(*_: object, **__: object) -> None:
        raise _StubError("Input CSV is missing required columns: foo")

    monkeypatch.setattr(get_target_data, "_IUPHAR_IMPORT_ERROR", None)
    monkeypatch.setattr(
        get_target_data,
        "iuphar_pp",
        SimpleNamespace(
            process_iuphar_targets=_raise_error,
            IUPHARPostProcessingError=_StubError,
        ),
    )

    result = get_target_data._postprocess_iuphar_export(source)

    assert result is None
    assert (
        "warning",
        "target_iuphar_postprocess_missing_columns",
        {"path": str(source), "error": "Input CSV is missing required columns: foo"},
    ) in logger_stub.events
