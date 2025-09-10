import logging
import sys
from pathlib import Path

import pytest

import get_activity_data as gad
import get_assay_data as gas
import get_document_data as gdd
import get_input_initialisation as gii
import get_target_data as gtd
import get_testitem_data as gti
import mapper_main as mp
import table_quality_main as tq
import get_document_type as gdt

CLI_MODULES = [
    (gad, []),
    (gas, []),
    (gdd, ["pubmed"]),
    (gii, []),
    (gtd, ["uniprot"]),
    (gti, []),
    (mp, []),
    (tq, ["--table-name", "demo"]),
]


@pytest.mark.parametrize("module,args", CLI_MODULES)
def test_main_invalid_config_returns_nonzero(tmp_path: Path, module, args, caplog):
    """Ensure CLI exits with non-zero code for malformed configuration."""
    bad_cfg = tmp_path / "config.yaml"
    bad_cfg.write_text("jobs:\n  concurrency: bad\n", encoding="utf8")
    caplog.set_level(logging.ERROR)
    exit_code = module.main(["--config", str(bad_cfg)] + args)
    assert exit_code == 1
    assert "invalid configuration" in caplog.text


def test_document_type_invalid_config(tmp_path: Path, monkeypatch, caplog) -> None:
    """Special handling for get_document_type.main lacking argv parameter."""
    bad_cfg = tmp_path / "config.yaml"
    bad_cfg.write_text("jobs:\n  concurrency: bad\n", encoding="utf8")
    caplog.set_level(logging.ERROR)
    argv = [
        "prog",
        "--config",
        str(bad_cfg),
        "--input",
        "in.csv",
        "--output",
        "out.csv",
    ]
    monkeypatch.setattr(sys, "argv", argv)
    exit_code = gdt.main()
    assert exit_code == 1
    assert "invalid configuration" in caplog.text
