"""Unit tests for :mod:`scripts.get_target_data`."""

from __future__ import annotations

import argparse
import stat
from pathlib import Path
from typing import Callable, Iterable, Sequence

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

    def debug(self, event: str, **payload: object) -> None:
        self.events.append(("debug", event, dict(payload)))

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
    "filename, expected",
    [
        ("output.target_20251005.csv", True),
        ("targets_20251005_normalized.csv", True),
        ("output.targets_20251005.csv.tmp", True),
        (".output.targets_20251005.csv_normalized.tmp", True),
        ("output.targets_uniprot.csv", True),
        ("targets.csv", True),
        ("out.csv", False),
        ("out_chembl.csv", False),
        ("out_uniprot.csv", False),
        ("out_iuphar.csv", False),
        ("out_raw.csv", False),
        ("out_normalized.csv", False),
    ],
)
def test_is_supported_target_export__cases(
    filename: str, expected: bool, tmp_path: Path
) -> None:
    path = tmp_path / filename
    path.write_text("", encoding="utf-8")

    result = get_target_data._is_supported_target_export(path)

    assert result is expected


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


def test_prepare_raw_destination__removes_existing(tmp_path: Path) -> None:
    destination = tmp_path / "raw.csv"
    destination.write_text("old", encoding="utf-8")

    get_target_data._prepare_raw_destination(destination)

    assert not destination.exists()
    assert destination.parent.exists()


def test_prepare_raw_destination__handles_permission_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "raw.csv"
    destination.write_text("old", encoding="utf-8")

    unlink_calls = {"count": 0}
    original_unlink = get_target_data.Path.unlink

    def _fake_unlink(self: Path, *args: object, **kwargs: object) -> None:
        if self == destination and unlink_calls["count"] == 0:
            unlink_calls["count"] += 1
            raise PermissionError("denied")
        unlink_calls["count"] += 1
        original_unlink(self, *args, **kwargs)

    monkeypatch.setattr(
        get_target_data.Path, "unlink", _fake_unlink, raising=False
    )

    chmod_modes: list[int] = []

    def _fake_chmod(path: Path, mode: int) -> None:
        if path == destination:
            chmod_modes.append(mode)

    monkeypatch.setattr(get_target_data.os, "chmod", _fake_chmod)

    get_target_data._prepare_raw_destination(destination)

    assert unlink_calls["count"] >= 2
    expected_mode = stat.S_IRUSR | stat.S_IWUSR | getattr(stat, "S_IWRITE", 0)
    assert chmod_modes == [expected_mode]
    assert not destination.exists()


def test_prepare_raw_destination__raises_when_unlink_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "raw.csv"
    destination.write_text("old", encoding="utf-8")

    def _always_fail(*args: object, **kwargs: object) -> None:
        raise PermissionError("locked")

    monkeypatch.setattr(
        get_target_data.Path, "unlink", _always_fail, raising=False
    )
    monkeypatch.setattr(get_target_data.os, "chmod", lambda *_, **__: None)

    with pytest.raises(OSError):
        get_target_data._prepare_raw_destination(destination)


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

    args = argparse.Namespace(input_csv=input_csv, final_out=output_csv)

    exit_code = get_target_data.run_uniprot(cfg, args)

    assert exit_code == 0
    assert recorded["path"] == output_csv
    assert isinstance(recorded["context"], get_target_data.IsoformPostprocessContext)
    assert recorded["context"].args is args
    assert recorded["ambiguous"] is None


def test_run_uniprot__doc_quality_reports(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg.system.doc_quality.enable = True
    cfg.system.doc_quality.sample_rows = 5
    input_csv = tmp_path / "uniprot_ids.csv"
    input_csv.write_text("uniprot_id\nP12345\n", encoding="utf-8")
    output_csv = tmp_path / "output.target_20250101.csv"

    monkeypatch.setattr(get_target_data.uu, "init_session", lambda *_, **__: None)

    def _fake_process(**kwargs: object) -> None:
        Path(kwargs["output_csv"]).write_text(
            "uniprot_id\nP12345\n", encoding=cfg.io.csv_encoding
        )

    monkeypatch.setattr(get_target_data.uu, "process", _fake_process)

    monkeypatch.setattr(
        get_target_data,
        "_postprocess_target_exports",
        lambda *_, **__: None,
    )

    captured: dict[str, object] = {}

    def _fake_build_quality(
        quality_cfg: object,
        *,
        table_name: str | Path,
        destination: Path | str | None,
    ) -> Callable[[pd.DataFrame], None]:
        captured["sample_rows"] = getattr(quality_cfg, "sample_rows", None)
        captured["table_name"] = Path(table_name).name
        if destination is None:
            captured["destination_dir"] = None
        else:
            captured["destination_dir"] = Path(destination)

        def _hook(df: pd.DataFrame) -> None:
            captured["df"] = df.copy()
            captured["hook_called"] = True

        return _hook

    monkeypatch.setattr(
        get_target_data,
        "build_table_quality_hook",
        _fake_build_quality,
    )

    args = argparse.Namespace(input_csv=input_csv, final_out=output_csv)

    exit_code = get_target_data.run_uniprot(cfg, args)

    assert exit_code == 0
    assert captured["table_name"] == output_csv.resolve().stem
    assert captured["destination_dir"] == output_csv.resolve().parent
    assert captured.get("hook_called") is True
    pd.testing.assert_frame_equal(
        captured["df"], pd.DataFrame({"uniprot_id": ["P12345"]})
    )


def test_fetch_uniprot__no_candidates_writes_empty_output(
    cfg: Config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / "output_uniprot.csv"
    chembl_df = pd.DataFrame({"target_chembl_id": ["CHEMBL1"]})

    plan = get_target_data._UniprotQueryPlan(
        unique_records=[],
        row_candidates=[[] for _ in chembl_df.index],
        row_index=list(chembl_df.index),
        candidate_columns=[],
    )

    monkeypatch.setattr(
        get_target_data,
        "_build_uniprot_query_plan",
        lambda *_: plan,
    )

    def _unexpected(*_: object, **__: object) -> None:  # pragma: no cover - defensive
        raise AssertionError("run_uniprot should not be invoked when no candidates are found")

    monkeypatch.setattr(get_target_data, "run_uniprot", _unexpected)

    result = get_target_data.fetch_uniprot(cfg, chembl_df, output_csv)

    assert output_csv.exists()
    pd.testing.assert_frame_equal(
        result,
        pd.DataFrame(
            {
                "uniprot_id": pd.Series(dtype=object),
                "original_id": pd.Series(dtype=object),
                "source_column": pd.Series(dtype=object),
                "mapping_uniprot_id": pd.Series(dtype=object),
            }
        ),
        check_index_type=False,
    )

    written = pd.read_csv(
        output_csv, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
    )
    assert list(written.columns) == [
        "uniprot_id",
        "original_id",
        "source_column",
        "mapping_uniprot_id",
    ]
    assert written.empty


def test_fetch_uniprot__missing_output_creates_placeholder(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    logger_stub: _MemoryLogger,
) -> None:
    cfg.io.csv_sep = ","
    cfg.io.csv_encoding = "utf-8"
    cfg.target.all.uniprot_column = "uniprot_id"
    cfg.target.all.target_csv = tmp_path / "all_targets.csv"
    cfg.target.all.family_csv = tmp_path / "all_families.csv"
    cfg.target.iuphar.target_csv = tmp_path / "iuphar_targets.csv"
    cfg.target.iuphar.family_csv = tmp_path / "iuphar_families.csv"

    chembl_df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprot_id": ["P12345"],
            "mapping_uniprot_id": ["P12345"],
            "pref_name": ["Example"],
        }
    )
    uniprot_df = pd.DataFrame({"uniprot_id": ["P12345"], "pref_name": ["Example"]})

    output_csv = tmp_path / "output.targets_iuphar.csv"
    captured: dict[str, Path] = {}

    def _fake_run_iuphar(local_cfg: Config, args: argparse.Namespace) -> int:
        captured["cfg"] = local_cfg
        captured["output_csv"] = Path(args.output_csv)
        return 0

    monkeypatch.setattr(get_target_data, "run_iuphar", _fake_run_iuphar)

    combined_df, iuphar_df = get_target_data.fetch_iuphar(
        cfg, chembl_df, uniprot_df, output_csv
    )

    assert captured["cfg"] is cfg
    assert captured["output_csv"] == output_csv
    assert output_csv.exists()
    assert "uniprot_id" in combined_df.columns
    pd.testing.assert_frame_equal(
        iuphar_df,
        pd.DataFrame({"uniprot_id": pd.Series(dtype=object)}),
        check_index_type=False,
    )
    assert any(
        level == "warning"
        and event == "missing_iuphar_output_file"
        and payload["path"] == str(output_csv)
        for level, event, payload in logger_stub.events
    )

    output_csv = tmp_path / "output_uniprot.csv"
    chembl_df = pd.DataFrame(
        {"target_chembl_id": ["CHEMBL1"], "uniprot_id": ["P12345"]}
    )

    plan = get_target_data._UniprotQueryPlan(
        unique_records=[
            {
                "uniprot_id": "P12345",
                "original_id": "P12345",
                "source_column": "uniprot_id",
            }
        ],
        row_candidates=[
            [get_target_data._UniprotCandidate("P12345", "uniprot_id", "P12345")]
        ],
        row_index=list(chembl_df.index),
        candidate_columns=["uniprot_id"],
    )

    monkeypatch.setattr(
        get_target_data, "_build_uniprot_query_plan", lambda *_: plan
    )

    def _fake_run(cfg_obj: Config, args: argparse.Namespace) -> int:
        assert Path(args.final_out) == output_csv
        return 0

    monkeypatch.setattr(get_target_data, "run_uniprot", _fake_run)

    result = get_target_data.fetch_uniprot(cfg, chembl_df, output_csv)

    assert output_csv.exists()
    written = pd.read_csv(
        output_csv, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
    )
    assert list(written.columns) == [
        "uniprot_id",
        "original_id",
        "source_column",
        "mapping_uniprot_id",
    ]
    assert written.empty

    assert len(result.index) == 1
    assert result.loc[result.index[0], "uniprot_id"] == "P12345"
    assert any(
        level == "warning"
        and event == "fetch_uniprot_output_missing"
        and payload == {"path": str(output_csv)}
        for level, event, payload in logger_stub.events
    )


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


def test_postprocess_isoform_export__runs_for_custom_name(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "custom_export.csv"
    source.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    called: list[tuple[tuple[object, ...], dict[str, object]]] = []

    def _record(*args: object, **kwargs: object) -> str:
        called.append((args, dict(kwargs)))
        return str(source)

    monkeypatch.setattr(get_target_data.target_pp, "process_targets", _record)

    result = get_target_data._postprocess_isoform_export(source, cfg=cfg)

    assert result == source
    assert called == [((str(source),), {"verbose": True})]
    assert (
        "info",
        "target_isoform_postprocess_done",
        {"path": str(source), "source": str(source)},
    ) in logger_stub.events


@pytest.mark.parametrize(
    "filename",
    ["output.target_20250101.csv", "targets.csv"],
)
def test_postprocess_target_exports__chains_helpers(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    filename: str,
) -> None:
    source = tmp_path / filename
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

    def _fake_names(path: Path, *, cfg: Config) -> Path:
        assert Path(path) == source
        call_order.append("names")
        return source

    def _fake_iuphar(path: Path, *, verbose: bool = True) -> Path:
        assert Path(path) == source
        call_order.append("iuphar")
        return source

    monkeypatch.setattr(get_target_data, "_postprocess_organism_export", _fake_organism)
    monkeypatch.setattr(get_target_data, "_postprocess_isoform_export", _fake_isoform)
    monkeypatch.setattr(get_target_data, "_postprocess_names_export", _fake_names)
    monkeypatch.setattr(get_target_data, "_postprocess_iuphar_export", _fake_iuphar)

    get_target_data._postprocess_target_exports(source, cfg=cfg)

    assert call_order == ["organism", "isoform", "names", "iuphar"]


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
