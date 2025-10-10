"""CLI argument handling tests for get_target_data pipeline."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from library.cli import parser as cli_parser
from library.config import ConfigMetadata
from library.pipelines.target.defaults import TARGET_MODE_DEFAULTS
from scripts import get_target_data


class _ConfigStub(SimpleNamespace):
    pass


def _make_stub_config() -> _ConfigStub:
    def _mode_defaults(name: str) -> SimpleNamespace:
        defaults = TARGET_MODE_DEFAULTS[name]
        return SimpleNamespace(
            column=defaults.column,
            chunk_size=defaults.chunk_size,
            timeout=defaults.timeout,
            limit=defaults.limit,
            offset=defaults.offset,
        )

    io_cfg = SimpleNamespace(
        output_stamp_mode="omit",
        output_dir=Path("output").resolve(),
    )
    return _ConfigStub(
        target=SimpleNamespace(
            chembl=_mode_defaults("chembl"),
            uniprot=_mode_defaults("uniprot"),
            iuphar=_mode_defaults("iuphar"),
            all=_mode_defaults("all"),
        ),
        local=SimpleNamespace(io=io_cfg),
        io=io_cfg,
    )


@pytest.mark.unit
def test_resolve_target_parameters__chembl_cli_overrides(monkeypatch):
    parser, _ = get_target_data.build_parser()
    args = parser.parse_args(
        [
            "chembl",
            "--chunk-size",
            "7",
            "--timeout",
            "45",
            "--limit",
            "5",
            "--offset",
            "3",
        ]
    )
    cfg = _make_stub_config()

    entries = get_target_data._resolve_target_parameters("chembl", cfg, args)

    assert cfg.target.chembl.chunk_size == 7
    assert cfg.target.chembl.timeout == pytest.approx(45.0)
    assert cfg.target.chembl.limit == 5
    assert cfg.target.chembl.offset == 3
    cli_sources = {entry.name for entry in entries if entry.source == "cli"}
    assert cli_sources == {"chunk_size", "timeout", "limit", "offset"}
    default_column = [entry for entry in entries if entry.name == "column"]
    assert default_column and default_column[0].source == "default"


@pytest.mark.unit
def test_resolve_target_parameters__config_precedence():
    cfg = _make_stub_config()
    cfg.target.uniprot.chunk_size = 250
    cfg.target.uniprot.timeout = 60.0
    parser, _ = get_target_data.build_parser()
    args = parser.parse_args(["uniprot"])

    entries = get_target_data._resolve_target_parameters("uniprot", cfg, args)

    assert cfg.target.uniprot.chunk_size == 250
    assert cfg.target.uniprot.timeout == pytest.approx(60.0)
    sources = {entry.name: entry.source for entry in entries}
    assert sources["chunk_size"] == "config"
    assert sources["timeout"] == "config"


@pytest.mark.unit
def test_resolve_target_parameters__all_global_fallback():
    cfg = _make_stub_config()
    parser, _ = get_target_data.build_parser()
    args = parser.parse_args(
        [
            "all",
            "--chunk-size",
            "9",
            "--chembl-limit",
            "4",
        ]
    )

    get_target_data._resolve_target_parameters("all", cfg, args)

    assert cfg.target.all.chunk_size == 9
    assert cfg.target.chembl.chunk_size == 9
    assert cfg.target.chembl.limit == 4


@pytest.mark.unit
def test_prepare_io_paths__output_alias_sets_final_out(tmp_path):
    parser, _ = get_target_data.build_parser()
    args = parser.parse_args(
        [
            "all",
            "--base-path",
            str(tmp_path),
            "--output",
            "custom.csv",
        ]
    )

    get_target_data.prepare_io_paths(
        args,
        input_default=get_target_data.DEFAULT_INPUT_NAME,
        output_stem=get_target_data.DEFAULT_OUTPUT_STEM,
    )

    expected = (tmp_path / "custom.csv").resolve()
    assert args.final_out == expected
    assert args.output_csv == expected
    assert args.output_dir is None


@pytest.mark.unit
def test_prepare_io_paths__default_filename_without_date(tmp_path):
    parser, _ = get_target_data.build_parser()
    args = parser.parse_args(
        [
            "all",
            "--base-path",
            str(tmp_path),
        ]
    )

    get_target_data.prepare_io_paths(
        args,
        input_default=get_target_data.DEFAULT_INPUT_NAME,
        output_stem=get_target_data.DEFAULT_OUTPUT_STEM,
    )

    expected = (
        tmp_path / f"output.{get_target_data.DEFAULT_OUTPUT_STEM}.csv"
    ).resolve()
    assert args.final_out == expected
    assert getattr(args, "date", None) is None


@pytest.mark.unit
def test_prepare_io_paths__respects_explicit_date(tmp_path):
    parser, _ = get_target_data.build_parser()
    args = parser.parse_args(
        [
            "all",
            "--base-path",
            str(tmp_path),
            "--date",
            "20250228",
        ]
    )

    get_target_data.prepare_io_paths(
        args,
        input_default=get_target_data.DEFAULT_INPUT_NAME,
        output_stem=get_target_data.DEFAULT_OUTPUT_STEM,
    )

    expected = (
        tmp_path / f"output.{get_target_data.DEFAULT_OUTPUT_STEM}_20250228.csv"
    ).resolve()
    assert args.final_out == expected
    assert args.date == "20250228"


@pytest.mark.unit
def test_target_iuphar_defaults__align_with_cli() -> None:
    cfg = _make_stub_config()
    defaults = TARGET_MODE_DEFAULTS["iuphar"]

    assert cfg.target.iuphar.column == defaults.column
    assert cfg.target.iuphar.chunk_size == defaults.chunk_size
    assert cfg.target.iuphar.timeout == pytest.approx(defaults.timeout)
    assert cfg.target.iuphar.offset == defaults.offset


@pytest.mark.unit
def test_run_logs_parameter_sources__default_sources(monkeypatch, tmp_path):
    parser, _ = get_target_data.build_parser()
    args = parser.parse_args(["chembl"])
    cfg = _make_stub_config()

    input_path = tmp_path / "targets.csv"
    input_path.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")
    args.input_csv = input_path

    logged = []

    def fake_info(event: str, **kwargs: object) -> None:  # type: ignore[override]
        logged.append((event, kwargs))

    args.func = lambda *_: 0
    monkeypatch.setattr(get_target_data.logger, "info", fake_info)

    result = get_target_data.run(cfg, args)
    assert result == 0
    cli_logs = [kwargs for event, kwargs in logged if event == "cli_parameter"]
    assert cli_logs
    assert {entry["source"] for entry in cli_logs} == {"default"}


@pytest.mark.unit
def test_parser_validation__rejects_zero_chunk_size():
    parser, _ = get_target_data.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["chembl", "--chunk-size", "0"])


@pytest.mark.unit
def test_apply_config_overrides__missing_config_attribute(monkeypatch, tmp_path):
    parser, _ = get_target_data.build_parser()
    args = parser.parse_args(["all"])
    subparser = parser.subparsers_map["all"]
    dummy_config = tmp_path / "config.yaml"
    dummy_config.write_text("{}", encoding="utf-8")

    def fake_load_config(*_args, **_kwargs):
        cfg = _make_stub_config()
        cfg.target.iuphar.column = None
        metadata = ConfigMetadata(snapshot={}, sources={})
        return cfg, metadata

    monkeypatch.setattr(cli_parser, "load_config", fake_load_config)
    warnings: list[tuple[str, dict[str, object]]] = []

    def fake_warning(event: str, **kwargs: object) -> None:
        warnings.append((event, kwargs))

    monkeypatch.setattr(cli_parser.logger, "warning", fake_warning)

    cfg = cli_parser.apply_config_overrides(
        args,
        subparser,
        dummy_config,
        mapping={"iuphar_column": "target.iuphar.column"},
        base_parser=parser,
    )

    assert isinstance(cfg, _ConfigStub)
    assert args.iuphar_column == TARGET_MODE_DEFAULTS["iuphar"].column
    assert warnings
    matching = [
        payload
        for event, payload in warnings
        if payload.get("argument") == "iuphar_column"
    ]
    assert matching
    assert matching[0]["path"] == "sources.chembl.pipelines.target.iuphar.column"
    assert matching[0]["error"]
