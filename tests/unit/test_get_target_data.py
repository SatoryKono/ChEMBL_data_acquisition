"""Unit tests for :mod:`scripts.get_target_data`."""

from __future__ import annotations

import argparse
import stat
from collections.abc import Callable, Iterable
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from scripts import get_target_data

from library.config import Config
from library.postprocess.common import PostprocessingPipelineResult


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


def _events_of_level(
    logger: _MemoryLogger, level: str
) -> list[tuple[str, str, dict[str, object]]]:
    return [event for event in logger.events if event[0] == level]


@pytest.fixture()
def logger_stub(monkeypatch: pytest.MonkeyPatch) -> _MemoryLogger:
    logger = _MemoryLogger()
    monkeypatch.setattr(get_target_data, "logger", logger)
    return logger


@pytest.mark.parametrize(
    "latin, cyrillic",
    list(get_target_data._RUSSIAN_KEYBOARD_MAP.items())[:10],
)
def test_translate_keyboard_layout__single_characters(
    latin: str, cyrillic: str
) -> None:
    assert get_target_data._translate_keyboard_layout(latin) == cyrillic
    assert get_target_data._translate_keyboard_layout(latin.upper()) == cyrillic.upper()


@pytest.mark.parametrize("text", ["Hello", "Test-Value", "123", "aA"])
def test_translate_keyboard_layout__composite_strings(text: str) -> None:
    translated = get_target_data._translate_keyboard_layout(text)
    assert len(translated) == len(text)
    for original, result in zip(text, translated, strict=False):
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
    "filename, expected",
    [
        ("targets_normalized.csv", "targets.csv"),
    ],
)
def test_normalise_target_export_name__normalised_suffix_cases(
    filename: str, expected: str
) -> None:
    result = get_target_data._normalise_target_export_name(Path(filename))

    assert result == expected


def test_resolve_output_metadata__strips_working_prefix() -> None:
    table_name, date_tag = get_target_data._resolve_output_metadata(
        Path(".output.targets_20240101.csv.tmp")
    )

    assert table_name == "targets"
    assert date_tag == "20240101"


def test_resolve_output_metadata__normalizes_table_hint(tmp_path: Path) -> None:
    output_path = tmp_path / "custom.csv"

    table_name, date_tag = get_target_data._resolve_output_metadata(
        output_path,
        table_hint=".output.targets",
        date_hint="20251231",
    )

    assert table_name == "targets"
    assert date_tag == "20251231"


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


def test_prepare_raw_destination__removes_existing(tmp_path: Path, cfg: Config) -> None:
    destination = tmp_path / "raw.csv"
    destination.write_text("old", encoding="utf-8")

    get_target_data._prepare_raw_destination(destination, cfg=cfg)

    assert not destination.exists()
    assert destination.parent.exists()


def test_prepare_raw_destination__handles_permission_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
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

    monkeypatch.setattr(get_target_data.Path, "unlink", _fake_unlink, raising=False)

    chmod_modes: list[int] = []

    def _fake_chmod(path: Path, mode: int) -> None:
        if path == destination:
            chmod_modes.append(mode)

    monkeypatch.setattr(get_target_data.os, "chmod", _fake_chmod)

    get_target_data._prepare_raw_destination(destination, cfg=cfg)

    assert unlink_calls["count"] >= 2
    expected_mode = stat.S_IRUSR | stat.S_IWUSR | getattr(stat, "S_IWRITE", 0)
    assert chmod_modes == [expected_mode]
    assert not destination.exists()


@pytest.mark.unit
def test_run_chembl__updates_final_out_with_standard_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """The CLI runner should expose the dataset path generated by the pipeline."""

    final_out = tmp_path / ".output.targets_20240101.csv_chembl.csv"
    canonical_path = (
        tmp_path / "output..output.targets_20240101.csv_chembl_20240101.csv"
    )
    canonical_path.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    cfg.io.output_dir = tmp_path
    cfg.io.exist_ok = True

    monkeypatch.setattr(
        get_target_data.io,
        "read_ids",
        lambda *_, **__: iter(["CHEMBL1"]),
    )

    class _DummyWriter:
        def __init__(self, destination: Path, *, cfg: Config, reindex_columns: bool):
            self.destination = destination

        def write(self, frame: pd.DataFrame) -> None:  # pragma: no cover - stub
            return None

    monkeypatch.setattr(get_target_data, "_RawDumpStreamWriter", _DummyWriter)
    monkeypatch.setattr(
        get_target_data,
        "_finalize_raw_dump_writer",
        lambda *args, **kwargs: True,
    )

    from library.cli_utils import PipelineExecutionResult

    def _fake_pipeline(**_: object) -> PipelineExecutionResult:
        return PipelineExecutionResult(
            exit_code=0,
            dataset_path=canonical_path,
            failure_path=None,
            metadata_path=None,
        )

    monkeypatch.setattr(get_target_data, "_run_pipeline_with_meta", _fake_pipeline)

    args = argparse.Namespace(
        input_csv=tmp_path / "input.csv",
        final_out=final_out,
        output_csv=final_out,
        raw_out=None,
        raw_format="csv",
        id_cols=None,
        offset=0,
        normalize_at_export=True,
        no_reindex_raw=False,
    )

    exit_code = get_target_data.run_chembl(cfg, args)

    assert exit_code == 0
    assert args.final_out == canonical_path
    assert args.output_csv == canonical_path


@pytest.mark.unit
def test_fetch_chembl__uses_updated_final_out_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Fetch helper must read the dataset emitted by the pipeline."""

    final_out = tmp_path / ".output.targets_20240101.csv_chembl.csv"
    canonical_path = (
        tmp_path / "output..output.targets_20240101.csv_chembl_20240101.csv"
    )
    canonical_df = pd.DataFrame(
        [{"target_chembl_id": "CHEMBL1", "pref_name": "Example"}]
    )

    def _fake_run_chembl(cfg_arg: Config, args: argparse.Namespace) -> int:
        assert cfg_arg is cfg
        canonical_df.to_csv(
            canonical_path,
            index=False,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
        args.final_out = canonical_path
        args.output_csv = canonical_path
        return 0

    monkeypatch.setattr(get_target_data, "run_chembl", _fake_run_chembl)

    cfg.io.output_dir = tmp_path
    cfg.io.exist_ok = True

    df = get_target_data.fetch_chembl(
        cfg,
        input_csv=tmp_path / "input.csv",
        final_out=final_out,
    )

    assert df.equals(canonical_df)


@pytest.mark.unit
def test_fetch_chembl__cleans_up_standard_outputs_when_requested(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
    final_out = tmp_path / "output.targets_chembl.csv"

    cleanup_calls: list[Path] = []

    def _fake_run_chembl(cfg_obj: Config, args: argparse.Namespace) -> int:
        Path(args.final_out).write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
        return 0

    def _record_cleanup(path: Path) -> None:
        cleanup_calls.append(path)

    monkeypatch.setattr(get_target_data, "run_chembl", _fake_run_chembl)
    monkeypatch.setattr(
        get_target_data, "_cleanup_standard_output_artifacts", _record_cleanup
    )

    df = get_target_data.fetch_chembl(
        cfg,
        input_csv,
        final_out,
        cleanup_standard_outputs=True,
    )

    assert len(df.index) == 1
    assert cleanup_calls == [final_out]


@pytest.mark.unit
def test_cleanup_standard_output_artifacts__removes_expected_files(
    tmp_path: Path,
) -> None:
    base_output = tmp_path / "output.targets_chembl.csv"
    base_output.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    dataset = tmp_path / "output.targets_chembl_20240101.csv"
    dataset.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
    dataset_meta = Path(f"{dataset}.meta.yaml")
    dataset_meta.write_text("meta", encoding="utf-8")

    quality = tmp_path / "output.targets_chembl_20240101_quality_report_table.csv"
    quality.write_text("quality", encoding="utf-8")
    quality_meta = Path(f"{quality}.meta.yaml")
    quality_meta.write_text("meta", encoding="utf-8")

    corr = tmp_path / "output.targets_chembl_20240101_data_correlation_report_table.csv"
    corr.write_text("corr", encoding="utf-8")
    corr_meta = Path(f"{corr}.meta.yaml")
    corr_meta.write_text("meta", encoding="utf-8")

    get_target_data._cleanup_standard_output_artifacts(base_output)

    assert base_output.exists()
    assert not dataset.exists()
    assert not dataset_meta.exists()
    assert not quality.exists()
    assert not quality_meta.exists()
    assert not corr.exists()
    assert not corr_meta.exists()


def test_prepare_raw_destination__raises_when_unlink_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    destination = tmp_path / "raw.csv"
    destination.write_text("old", encoding="utf-8")

    def _always_fail(*args: object, **kwargs: object) -> None:
        raise PermissionError("locked")

    monkeypatch.setattr(get_target_data.Path, "unlink", _always_fail, raising=False)
    monkeypatch.setattr(get_target_data.os, "chmod", lambda *_, **__: None)

    with pytest.raises(OSError):
        get_target_data._prepare_raw_destination(destination, cfg=cfg)


def test_prepare_raw_destination__fails_when_parent_missing_and_exist_ok_false(
    tmp_path: Path, cfg: Config
) -> None:
    cfg.io.exist_ok = False
    destination = tmp_path / "missing" / "raw.csv"

    with pytest.raises(FileNotFoundError):
        get_target_data._prepare_raw_destination(destination, cfg=cfg)


def test_restore_legacy_output__restores_missing_file(
    tmp_path: Path, cfg: Config, logger_stub: _MemoryLogger
) -> None:
    source = tmp_path / "output" / "canonical.csv"
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_text("value\n", encoding="utf-8")

    destination = tmp_path / "final" / "legacy.csv"

    get_target_data._restore_legacy_output(destination, source, cfg=cfg)

    assert destination.exists()
    assert destination.read_text(encoding="utf-8") == source.read_text(encoding="utf-8")
    restored_events = [
        payload
        for level, event, payload in logger_stub.events
        if level == "warning" and event == "legacy_output_restored"
    ]
    assert restored_events
    assert restored_events[0]["source"] == str(source)
    assert restored_events[0]["destination"] == str(destination)


def test_restore_legacy_output__skips_when_destination_exists(
    tmp_path: Path, cfg: Config, logger_stub: _MemoryLogger
) -> None:
    source = tmp_path / "output.csv"
    source.write_text("new\n", encoding="utf-8")

    destination = tmp_path / "output.csv"
    destination.write_text("existing\n", encoding="utf-8")

    get_target_data._restore_legacy_output(destination, source, cfg=cfg)

    assert destination.read_text(encoding="utf-8") == "existing\n"
    restored_events = [
        event for _, event, _ in logger_stub.events if event == "legacy_output_restored"
    ]
    assert not restored_events


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


def test_ensure_merge_column_present__uses_alias(
    cfg: Config, logger_stub: _MemoryLogger
) -> None:
    cfg.target.all.uniprot_column = "canonical_uniprot"
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "mapping_uniprot_id": ["P12345"],
        }
    )

    result = get_target_data._ensure_merge_column_present(
        frame, cfg.target.all.uniprot_column, cfg
    )

    assert "canonical_uniprot" in result.columns
    assert result.loc[0, "canonical_uniprot"] == "P12345"
    assert any(
        event == "uniprot_merge_column_alias" for _, event, _ in logger_stub.events
    )


def test_ensure_merge_column_present__raises(
    cfg: Config, logger_stub: _MemoryLogger
) -> None:
    cfg.target.all.uniprot_column = "canonical_uniprot"
    frame = pd.DataFrame({"target_chembl_id": ["CHEMBL1"]})

    with pytest.raises(get_target_data.PipelineError):
        get_target_data._ensure_merge_column_present(
            frame, cfg.target.all.uniprot_column, cfg
        )

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
    assert (
        "info",
        "pipeline_skip_existing",
        {"output": str(final_out)},
    ) in logger_stub.events


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


def test_fetch_chembl__falls_back_to_normalized_output(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    logger_stub: _MemoryLogger,
) -> None:
    input_csv = tmp_path / "identifiers.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
    expected_output = tmp_path / "targets_chembl.csv"
    normalized_output = expected_output.with_name(
        f"{expected_output.stem}_normalized{expected_output.suffix}"
    )
    normalized_output.write_text(
        "target_chembl_id\nCHEMBL1\n", encoding=cfg.io.csv_encoding
    )

    def _fake_run(cfg_arg: Config, args: argparse.Namespace) -> int:
        assert args.final_out == expected_output
        return 0

    monkeypatch.setattr(get_target_data, "run_chembl", _fake_run)

    frame = get_target_data.fetch_chembl(cfg, input_csv, expected_output)

    assert list(frame["target_chembl_id"]) == ["CHEMBL1"]
    assert (
        "warning",
        "fetch_chembl_missing_expected_output",
        {
            "expected": str(expected_output),
            "fallback": str(normalized_output),
        },
    ) in logger_stub.events


def test_fetch_chembl__raises_when_output_missing(
    cfg: Config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "identifiers.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
    expected_output = tmp_path / "targets_chembl.csv"

    def _fake_run(cfg_arg: Config, args: argparse.Namespace) -> int:
        assert args.final_out == expected_output
        return 0

    monkeypatch.setattr(get_target_data, "run_chembl", _fake_run)

    with pytest.raises(FileNotFoundError):
        get_target_data.fetch_chembl(cfg, input_csv, expected_output)


def test_run_chembl__restores_final_output_without_legacy(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    logger_stub: _MemoryLogger,
) -> None:
    input_csv = tmp_path / "identifiers.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
    final_out = tmp_path / "output.targets_chembl.csv"
    dataset_dir = tmp_path / "standard"
    dataset_dir.mkdir()
    dataset_path = dataset_dir / "output.targets_chembl.csv"
    dataset_content = "target_chembl_id\nCHEMBL1\n"
    dataset_path.write_text(dataset_content, encoding=cfg.io.csv_encoding)

    monkeypatch.setattr(
        get_target_data.io,
        "read_ids",
        lambda path, column, cfg: iter(["CHEMBL1"]),
    )

    class _StubWriter:
        def write(self, chunk: pd.DataFrame) -> None:  # pragma: no cover - stub
            return None

        def finalize(self) -> None:
            return None

    monkeypatch.setattr(
        get_target_data,
        "_RawDumpStreamWriter",
        lambda *_, **__: _StubWriter(),
    )

    def _fake_run_pipeline(**kwargs: object) -> SimpleNamespace:
        assert Path(kwargs["output_path"]) == final_out
        return SimpleNamespace(exit_code=0, dataset_path=str(dataset_path))

    monkeypatch.setattr(
        get_target_data,
        "_run_pipeline_with_meta",
        _fake_run_pipeline,
    )

    monkeypatch.setattr(
        get_target_data,
        "build_table_quality_hook",
        lambda *_, **__: lambda *_: None,
    )

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=final_out,
        raw_out=None,
        raw_format="csv",
        no_reindex_raw=False,
        emit_legacy_artifacts=False,
    )

    exit_code = get_target_data.run_chembl(cfg, args)

    assert exit_code == 0
    assert final_out.read_text(encoding=cfg.io.csv_encoding) == dataset_content
    assert (
        "info",
        "chembl_final_output_restored",
        {
            "source": str(dataset_path),
            "destination": str(final_out),
        },
    ) in logger_stub.events


@pytest.mark.unit
def test_run_chembl__forces_standard_outputs_when_outputs_disabled(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    logger_stub: _MemoryLogger,
) -> None:
    input_csv = tmp_path / "identifiers.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
    final_out = tmp_path / "output.targets_chembl.csv"
    dataset_dir = tmp_path / "standard"
    dataset_dir.mkdir()
    dataset_path = dataset_dir / "output.targets_chembl.csv"
    dataset_path.write_text("target_chembl_id\nCHEMBL1\n", encoding=cfg.io.csv_encoding)

    def _fake_read_ids(path: Path, *, column: str, **_: object) -> Iterable[str]:
        assert column == cfg.target.chembl.column
        return iter(["CHEMBL1"])

    monkeypatch.setattr(get_target_data.io, "read_ids", _fake_read_ids)

    monkeypatch.setattr(
        get_target_data,
        "build_table_quality_hook",
        lambda *_, **__: lambda *_: None,
    )

    captured: dict[str, object] = {}

    def _fake_run_pipeline(**kwargs: object) -> SimpleNamespace:
        captured["emit_standard_outputs"] = kwargs["emit_standard_outputs"]
        captured["emit_legacy_artifacts"] = kwargs["emit_legacy_artifacts"]
        return SimpleNamespace(exit_code=0, dataset_path=str(dataset_path))

    monkeypatch.setattr(
        get_target_data,
        "_run_pipeline_with_meta",
        _fake_run_pipeline,
    )

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=final_out,
        output_csv=final_out,
        raw_out=None,
        raw_format="csv",
        id_cols=None,
        offset=0,
        normalize_at_export=True,
        no_reindex_raw=False,
        emit_standard_outputs=False,
        emit_legacy_artifacts=False,
    )

    exit_code = get_target_data.run_chembl(cfg, args)

    assert exit_code == 0
    assert captured["emit_standard_outputs"] is True
    assert captured["emit_legacy_artifacts"] is False
    assert (
        "info",
        "chembl_pipeline_outputs_forced",
        {"reason": "run_pipeline_requires_outputs"},
    ) in logger_stub.events
    assert args.emit_standard_outputs is True
    assert args.emit_legacy_artifacts is False


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
        raise AssertionError(
            "run_uniprot should not be invoked when no candidates are found"
        )

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

    monkeypatch.setattr(get_target_data, "_build_uniprot_query_plan", lambda *_: plan)

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


@pytest.mark.unit
def test_run_all__enables_standard_outputs_for_chembl(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
    final_out = tmp_path / "output.targets.csv"

    captured: dict[str, object] = {}

    def _fake_fetch_chembl(
        cfg_obj: Config,
        input_path: Path,
        final_path: Path,
        **kwargs: object,
    ) -> pd.DataFrame:
        captured["emit_standard_outputs"] = kwargs["emit_standard_outputs"]
        captured["cleanup_standard_outputs"] = kwargs["cleanup_standard_outputs"]
        return pd.DataFrame({"target_chembl_id": ["CHEMBL1"]})

    def _fake_fetch_uniprot(
        cfg_obj: Config,
        chembl_df: pd.DataFrame,
        output_path: Path,
    ) -> pd.DataFrame:
        return chembl_df.assign(uniprot_id=["P12345"])

    def _fake_fetch_iuphar(
        cfg_obj: Config,
        chembl_df: pd.DataFrame,
        uniprot_df: pd.DataFrame,
        output_path: Path,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        return uniprot_df, pd.DataFrame()

    def _fake_validate_and_write(
        merged: pd.DataFrame,
        destination: Path,
        cfg_obj: Config,
        **_: object,
    ) -> int:
        destination.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(get_target_data, "fetch_chembl", _fake_fetch_chembl)
    monkeypatch.setattr(get_target_data, "fetch_uniprot", _fake_fetch_uniprot)
    monkeypatch.setattr(get_target_data, "fetch_iuphar", _fake_fetch_iuphar)
    monkeypatch.setattr(get_target_data, "merge_results", lambda a, *_: a)
    monkeypatch.setattr(get_target_data, "validate_and_write", _fake_validate_and_write)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=final_out,
        raw_out=None,
        raw_format="csv",
        no_reindex_raw=False,
        disable_gtop=False,
        emit_legacy_artifacts=False,
        chembl_out=None,
        uniprot_out=None,
        iuphar_out=None,
        id_cols=None,
        date=None,
        _auto_output_stem=None,
    )

    exit_code = get_target_data.run_all(cfg, args)

    assert exit_code == 0
    assert captured.get("emit_standard_outputs") is True
    assert captured.get("cleanup_standard_outputs") is True
    assert final_out.exists()


@pytest.mark.unit
def test_run_chembl__coerces_none_emit_flag(
    cfg: Config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``None`` should be treated as ``True`` for ``emit_standard_outputs``."""

    final_out = tmp_path / "output.targets_chembl.csv"

    monkeypatch.setattr(
        get_target_data.io,
        "read_ids",
        lambda *_, **__: iter(["CHEMBL1"]),
    )

    class _StubWriter:
        def __init__(
            self, destination: Path, *, cfg: Config, reindex_columns: bool
        ) -> None:
            self.destination = destination

        def write(self, chunk: pd.DataFrame) -> None:  # pragma: no cover - stub
            return None

        def finalize(self) -> None:  # pragma: no cover - stub
            return None

    monkeypatch.setattr(get_target_data, "_RawDumpStreamWriter", _StubWriter)
    monkeypatch.setattr(
        get_target_data,
        "_finalize_raw_dump_writer",
        lambda *_, **__: True,
    )

    from library.cli_utils import PipelineExecutionResult

    captured: dict[str, object] = {}

    def _fake_pipeline(**kwargs: object) -> PipelineExecutionResult:
        captured["emit_standard_outputs"] = kwargs["emit_standard_outputs"]
        captured["emit_legacy_artifacts"] = kwargs["emit_legacy_artifacts"]
        final_out.write_text(
            "target_chembl_id\nCHEMBL1\n", encoding=cfg.io.csv_encoding
        )
        return PipelineExecutionResult(
            exit_code=0,
            dataset_path=final_out,
            failure_path=None,
            metadata_path=None,
        )

    monkeypatch.setattr(get_target_data, "_run_pipeline_with_meta", _fake_pipeline)

    args = argparse.Namespace(
        input_csv=tmp_path / "input.csv",
        final_out=final_out,
        output_csv=final_out,
        raw_out=None,
        raw_format="csv",
        id_cols=None,
        no_reindex_raw=False,
        normalize_at_export=True,
        emit_standard_outputs=None,
        date=None,
    )

    exit_code = get_target_data.run_chembl(cfg, args)

    assert exit_code == 0
    assert captured.get("emit_standard_outputs") is True
    assert captured.get("emit_legacy_artifacts") is False
    assert final_out.exists()


@pytest.mark.unit
def test_run_all__cleans_up_default_intermediate_outputs(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
    final_out = tmp_path / "output.targets.csv"

    created_paths: dict[str, Path] = {}

    def _fake_fetch_chembl(
        cfg_obj: Config,
        input_path: Path,
        final_path: Path,
        **kwargs: object,
    ) -> pd.DataFrame:
        final_path.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
        created_paths["chembl"] = final_path
        return pd.DataFrame({"target_chembl_id": ["CHEMBL1"]})

    def _fake_fetch_uniprot(
        cfg_obj: Config,
        chembl_df: pd.DataFrame,
        output_path: Path,
    ) -> pd.DataFrame:
        output_path.write_text(
            "target_chembl_id,uniprot_id\nCHEMBL1,P12345\n",
            encoding="utf-8",
        )
        created_paths["uniprot"] = output_path
        return chembl_df.assign(uniprot_id=["P12345"])

    def _fake_fetch_iuphar(
        cfg_obj: Config,
        chembl_df: pd.DataFrame,
        uniprot_df: pd.DataFrame,
        output_path: Path,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        output_path.write_text(
            "target_chembl_id,iuphar_id\nCHEMBL1,GTOP\n",
            encoding="utf-8",
        )
        created_paths["iuphar"] = output_path
        return uniprot_df, pd.DataFrame({"target_chembl_id": ["CHEMBL1"]})

    def _fake_merge_results(
        combined: pd.DataFrame, iuphar_df: pd.DataFrame, cfg_obj: Config
    ) -> pd.DataFrame:
        return combined.assign(iuphar_id=["GTOP"])

    def _fake_validate_and_write(
        merged: pd.DataFrame,
        destination: Path,
        cfg_obj: Config,
        **_: object,
    ) -> int:
        destination.write_text(
            "target_chembl_id,uniprot_id,iuphar_id\nCHEMBL1,P12345,GTOP\n",
            encoding="utf-8",
        )
        return 0

    monkeypatch.setattr(get_target_data, "fetch_chembl", _fake_fetch_chembl)
    monkeypatch.setattr(get_target_data, "fetch_uniprot", _fake_fetch_uniprot)
    monkeypatch.setattr(get_target_data, "fetch_iuphar", _fake_fetch_iuphar)
    monkeypatch.setattr(get_target_data, "merge_results", _fake_merge_results)
    monkeypatch.setattr(get_target_data, "validate_and_write", _fake_validate_and_write)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=final_out,
        raw_out=None,
        raw_format="csv",
        no_reindex_raw=False,
        disable_gtop=False,
        emit_legacy_artifacts=False,
        chembl_out=None,
        uniprot_out=None,
        iuphar_out=None,
        id_cols=None,
        date=None,
        _auto_output_stem=None,
    )

    exit_code = get_target_data.run_all(cfg, args)

    chembl_path = final_out.with_name(final_out.stem + "_chembl.csv")
    uniprot_path = final_out.with_name(final_out.stem + "_uniprot.csv")
    iuphar_path = final_out.with_name(final_out.stem + "_iuphar.csv")

    assert exit_code == 0
    assert final_out.exists()
    assert created_paths["chembl"] == chembl_path
    assert created_paths["uniprot"] == uniprot_path
    assert created_paths["iuphar"] == iuphar_path
    assert not chembl_path.exists()
    assert not uniprot_path.exists()
    assert not iuphar_path.exists()


def test_run_target_postprocess_if_requested__disabled_flag_logs_skip(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
) -> None:
    source = tmp_path / "output.targets_20250101.csv"
    source.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    result = get_target_data.run_target_postprocess_if_requested(
        source,
        cfg=cfg,
        args=argparse.Namespace(postprocess=False),
    )

    assert result is None
    assert (
        "info",
        "Postprocessing skipped (flag --postprocess not set)",
        {},
    ) in logger_stub.events


def test_run_target_postprocess_if_requested__skips_unsupported_export(
    cfg: Config,
    tmp_path: Path,
    logger_stub: _MemoryLogger,
) -> None:
    source = tmp_path / "output.targets_20250101_uniprot.csv"
    source.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    result = get_target_data.run_target_postprocess_if_requested(
        source,
        cfg=cfg,
        args=argparse.Namespace(postprocess=True),
    )

    assert result is None
    assert (
        "info",
        "target_postprocess_skipped",
        {"path": str(source), "reason": "unsupported_export_name"},
    ) in logger_stub.events


def test_run_target_postprocess_if_requested__invokes_pipeline(
    cfg: Config,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    logger_stub: _MemoryLogger,
) -> None:
    import library.postprocess as postprocess_module
    import library.postprocessing.targets as targets_postprocess

    source = tmp_path / "output.targets_20250101.csv"
    destination = source.with_name("output_postprocessed.targets.csv")
    source.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    report_path = destination.parent / "targets.postprocess.report.json"

    class _Result:
        def __init__(self) -> None:
            self.output_path = destination
            self.report_path = report_path
            self.metrics = SimpleNamespace(
                summary=lambda: {
                    "rows": 1,
                    "columns": 2,
                    "steps": 3,
                    "duration_s": 0.5,
                },
                pipeline_version="1.0.0",
                validation=None,
            )

    pipeline_calls: list[tuple[tuple[object, ...], dict[str, object]]] = []

    def _fake_run_postprocessing_pipeline(*args: object, **kwargs: object) -> _Result:
        pipeline_calls.append((args, dict(kwargs)))
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            "target_chembl_id,pref_name\nCHEMBL1,Alpha\n",
            encoding="utf-8",
        )
        report_path.write_text("{}", encoding="utf-8")
        return _Result()

    class _RuntimeConfig:
        def __init__(self, **kwargs: object) -> None:
            self.kwargs = kwargs

    monkeypatch.setattr(
        postprocess_module,
        "PostprocessingPipelineConfig",
        _RuntimeConfig,
    )
    monkeypatch.setattr(
        postprocess_module,
        "get_pipeline_config",
        lambda table, override: SimpleNamespace(params={}, pipeline_version="1.0.0"),
    )
    monkeypatch.setattr(
        postprocess_module,
        "get_csv_runtime_config",
        lambda *_: SimpleNamespace(separator=",", encoding="utf-8", chunksize=10000),
    )
    monkeypatch.setattr(
        postprocess_module,
        "run_postprocessing_pipeline",
        _fake_run_postprocessing_pipeline,
    )

    def _dummy_runner(
        df: object, *, pipeline_version: object | None = None, logger=None
    ):
        return df, SimpleNamespace(
            summary=lambda: {
                "rows": 1,
                "columns": 2,
                "steps": 3,
                "duration_s": 0.5,
            },
            pipeline_version=pipeline_version or "1.0.0",
            validation=None,
        )

    monkeypatch.setattr(targets_postprocess, "run_target_pipeline", _dummy_runner)
    monkeypatch.setattr(targets_postprocess, "validate_targets", lambda df, **_: df)

    result = get_target_data.run_target_postprocess_if_requested(
        source,
        cfg=cfg,
        args=argparse.Namespace(postprocess=True),
    )

    assert isinstance(result, _Result)
    assert pipeline_calls and pipeline_calls[0][0][0] == "targets"
    assert pipeline_calls[0][0][1] == source
    assert pipeline_calls[0][0][2] == destination

    events = [
        event for event in logger_stub.events if event[1] == "target_postprocess_done"
    ]
    assert events
    _, _, payload = events[0]
    assert payload["path"] == str(destination)
    assert payload["postprocess_rows"] == 1
    assert payload["postprocess_columns"] == 2
    assert payload["postprocess_steps"] == 3
    assert payload["postprocess_duration_s"] == 0.5


def test_finalize_raw_dump_writer__missing_method(
    tmp_path: Path, logger_stub: _MemoryLogger
) -> None:
    class _Writer:
        pass

    writer = _Writer()
    destination = tmp_path / "raw.csv"

    result = get_target_data._finalize_raw_dump_writer(
        writer,
        logger=logger_stub,
        destination=destination,
    )

    assert result is True
    debug_events = _events_of_level(logger_stub, "debug")
    assert (
        "debug",
        "raw_dump_finalize_missing",
        {"writer_type": "_Writer"},
    ) in debug_events


def test_finalize_raw_dump_writer__error(
    tmp_path: Path, logger_stub: _MemoryLogger
) -> None:
    class _Writer:
        def finalize(self) -> None:
            raise OSError("boom")

    writer = _Writer()
    destination = tmp_path / "raw.csv"

    result = get_target_data._finalize_raw_dump_writer(
        writer,
        logger=logger_stub,
        destination=destination,
    )

    assert result is False
    error_events = _events_of_level(logger_stub, "error")
    assert len(error_events) == 1
    level, event, payload = error_events[0]
    assert level == "error"
    assert event == "raw_dump_failed"
    assert payload["error"] == "boom"
    assert payload["path"] == str(destination)
    assert isinstance(payload["exc_info"], OSError)


def test_finalize_raw_dump_writer__success(
    tmp_path: Path, logger_stub: _MemoryLogger
) -> None:
    calls: list[str] = []

    class _Writer:
        def finalize(self) -> None:
            calls.append("finalize")

    writer = _Writer()
    destination = tmp_path / "raw.csv"

    result = get_target_data._finalize_raw_dump_writer(
        writer,
        logger=logger_stub,
        destination=destination,
    )

    assert result is True
    assert calls == ["finalize"]
    assert not _events_of_level(logger_stub, "error")


def test_raw_dump_stream_writer_finalize__creates_empty_csv(
    tmp_path: Path, cfg: Config
) -> None:
    destination = tmp_path / "raw.csv"
    writer = get_target_data._RawDumpStreamWriter(
        destination, cfg=cfg, reindex_columns=False
    )

    result = writer.finalize()

    assert result == destination
    assert destination.exists()
    content = destination.read_text(encoding="utf-8")
    assert content.lstrip("\ufeff").replace("\r", "") == "\n"


def test_raw_dump_stream_writer_write__permission_error_retries(
    tmp_path: Path,
    cfg: Config,
    monkeypatch: pytest.MonkeyPatch,
    logger_stub: _MemoryLogger,
) -> None:
    destination = tmp_path / "raw.csv"
    writer = get_target_data._RawDumpStreamWriter(
        destination, cfg=cfg, reindex_columns=False
    )

    frame = pd.DataFrame({"col": [1]})
    calls: list[str] = []
    original_to_csv = pd.DataFrame.to_csv

    def _flaky_to_csv(self, *args, **kwargs):  # type: ignore[override]
        calls.append("call")
        if len(calls) == 1:
            raise PermissionError("file is locked")
        return original_to_csv(self, *args, **kwargs)

    sleeps: list[float] = []

    monkeypatch.setattr(pd.DataFrame, "to_csv", _flaky_to_csv)
    monkeypatch.setattr(get_target_data.time, "sleep", sleeps.append)

    writer.write(frame)

    assert calls == ["call", "call"]
    assert sleeps == [
        get_target_data._RawDumpStreamWriter._CSV_WRITE_RETRY_BASE_SLEEP_S
    ]
    assert (
        destination.read_text(encoding="utf-8").lstrip("\ufeff").replace("\r", "")
        == "col\n1\n"
    )
    warning_events = _events_of_level(logger_stub, "warning")
    assert (
        "warning",
        "raw_dump_write_retry",
        {
            "attempt": 1,
            "wait_seconds": get_target_data._RawDumpStreamWriter._CSV_WRITE_RETRY_BASE_SLEEP_S,
            "path": str(destination),
            "error": "file is locked",
        },
    ) in warning_events


def test_validate_and_write__omits_failure_artifacts_when_legacy_disabled(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        get_target_data,
        "generate_qc_report",
        lambda *_, **__: pd.DataFrame(),
    )
    monkeypatch.setattr(
        get_target_data,
        "generate_correlation_report",
        lambda *_, **__: pd.DataFrame(),
    )

    cfg.io.output_dir = tmp_path
    output_path = tmp_path / "output.targets_20230101.csv"
    failure_path = tmp_path / "output.targets_20230101_failure_cases.csv"
    failure_meta = failure_path.with_name(failure_path.name + ".meta.yaml")
    failure_path.write_text("stale", encoding="utf-8")
    failure_meta.write_text("metadata", encoding="utf-8")

    df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "sequence_length": [123],
        }
    )

    exit_code = get_target_data.validate_and_write(
        df,
        output_path,
        cfg,
        emit_legacy_artifacts=False,
    )

    assert exit_code == 1
    assert output_path.exists()
    assert not failure_path.exists()
    assert not failure_meta.exists()


def test_validate_and_write__removes_postprocess_sidecars(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg.io.output_dir = tmp_path
    output_path = tmp_path / "output.targets_20230101.csv"

    source_df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "protein_class_pred_rule_id": ["ambiguous"],
        }
    )

    processed_df = source_df.assign(protein_class_pred_rule_id=["resolved"])

    monkeypatch.setattr(get_target_data, "normalize_targets", lambda df: df)
    monkeypatch.setattr(get_target_data, "add_pipeline_metadata", lambda df: df)
    monkeypatch.setattr(
        get_target_data,
        "_prepare_targets_for_schema",
        lambda df: (df, set(), set()),
    )

    class _ValidationStub:
        def __init__(self, data: pd.DataFrame) -> None:
            self.data = data
            self.failure_cases = pd.DataFrame()

    monkeypatch.setattr(
        get_target_data,
        "validate_targets",
        lambda df, return_result=True: _ValidationStub(df),
    )
    monkeypatch.setattr(
        get_target_data,
        "generate_qc_report",
        lambda *_args, **_kwargs: pd.DataFrame({"metric": [1]}),
    )
    monkeypatch.setattr(
        get_target_data,
        "generate_correlation_report",
        lambda *_args, **_kwargs: pd.DataFrame({"metric": [1]}),
    )
    monkeypatch.setattr(
        get_target_data,
        "build_table_quality_hook",
        lambda *args, **kwargs: lambda _df: None,
    )

    postprocess_output = tmp_path / "output_postprocessed.targets.csv"
    postprocess_report = tmp_path / "targets.postprocess.report.json"

    def _fake_postprocess(
        *args: object, **kwargs: object
    ) -> PostprocessingPipelineResult:
        postprocess_output.write_text("postprocessed", encoding="utf-8")
        postprocess_report.write_text("{}", encoding="utf-8")
        return PostprocessingPipelineResult(
            dataframe=processed_df,
            metrics=None,
            output_path=postprocess_output,
            report_path=postprocess_report,
        )

    monkeypatch.setattr(
        get_target_data,
        "run_target_postprocess_if_requested",
        _fake_postprocess,
    )

    context = get_target_data.IsoformPostprocessContext(
        args=argparse.Namespace(postprocess=True)
    )

    exit_code = get_target_data.validate_and_write(
        source_df,
        output_path,
        cfg,
        postprocess_context=context,
    )

    assert exit_code == 0
    assert not postprocess_output.exists()
    assert not postprocess_report.exists()

    dataset = pd.read_csv(output_path)
    pd.testing.assert_frame_equal(
        dataset[processed_df.columns].reset_index(drop=True),
        processed_df.reset_index(drop=True),
    )

    produced = sorted(path.name for path in tmp_path.iterdir())
    assert produced == [
        "output.targets_20230101.csv",
        "output.targets_20230101.meta.yaml",
        "output.targets_20230101_data_correlation_report_table.csv",
        "output.targets_20230101_quality_report_table.csv",
    ]
