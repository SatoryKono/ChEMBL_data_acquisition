from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from library.cli.commands import get_data


def _make_run_config(tmp_path: Path) -> get_data.PipelineRunConfig:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    config_path = base_path / "config.yaml"
    config_path.write_text("io: {}\n", encoding="utf-8")
    input_files = get_data.PipelineInputFiles.from_mapping({"document": "document.csv"})
    output_stems = get_data.PipelineOutputStems.from_mapping({"document": "document"})
    subcommands = get_data.PipelineSubcommands.from_mapping({"document": None})
    return get_data.PipelineRunConfig(
        base_path=base_path,
        input_dir=input_dir,
        output_dir=output_dir,
        config_path=config_path,
        date_prefix="20200101",
        log_level="INFO",
        limit=None,
        force=False,
        skip_existing=False,
        dry_run=False,
        input_files=input_files,
        output_stems=output_stems,
        subcommands=subcommands,
    )


@pytest.mark.unit
def test_output_path__rejects_parent_escape(tmp_path: Path) -> None:
    cfg = _make_run_config(tmp_path)
    stems = cfg.output_stems.to_dict()
    stems["document"] = "../escape.csv"
    insecure_cfg = replace(
        cfg,
        output_stems=get_data.PipelineOutputStems.from_mapping(stems),
    )

    with pytest.raises(ValueError) as excinfo:
        insecure_cfg.output_path("document")

    message = str(excinfo.value)
    assert "escape" in message
    assert "output directory" in message
