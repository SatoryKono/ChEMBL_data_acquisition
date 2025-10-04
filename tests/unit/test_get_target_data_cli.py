"""CLI argument handling tests for get_target_data pipeline."""

from __future__ import annotations

import pytest

from library.config import Config
from scripts import get_target_data


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
    cfg = Config()

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
    cfg = Config()
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
    cfg = Config()
    parser, _ = get_target_data.build_parser()
    args = parser.parse_args([
        "all",
        "--chunk-size",
        "9",
        "--chembl-limit",
        "4",
    ])

    get_target_data._resolve_target_parameters("all", cfg, args)

    assert cfg.target.all.chunk_size == 9
    assert cfg.target.chembl.chunk_size == 9
    assert cfg.target.chembl.limit == 4


@pytest.mark.unit
def test_run_logs_parameter_sources__default_sources(monkeypatch, tmp_path):
    parser, _ = get_target_data.build_parser()
    args = parser.parse_args(["chembl"])
    cfg = Config()

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
