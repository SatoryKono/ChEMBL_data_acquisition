import io
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from pytest import LogCaptureFixture, MonkeyPatch

from library.cli import LoggerConfig, configure_logger
from library.logging_setup import Logger
from scripts import get_activity_data as gad
from scripts import get_assay_data as gas
from scripts import get_document_data as gdd
from scripts import get_document_type as gdoctype
from scripts import get_input_initialisation as gii
from scripts import get_target_data as gtd
from scripts import get_testitem_data as gtdt
from scripts import mapper_main as mapper
from scripts import table_quality_main as tqm

CLIS = [
    (gad.main, [], False),
    (gas.main, [], False),
    (gdd.main, ["chembl"], False),
    (gtdt.main, [], False),
    (mapper.main, [], False),
    (tqm.main, ["--table-name", "demo"], False),
    (gii.main, [], False),
    (gtd.main, ["chembl"], False),
    (gdoctype.main, ["--input", "in.csv", "--output", "out.csv"], True),
]


@pytest.mark.parametrize("entry, extra, use_sys", CLIS)
def test_malformed_config_exits(
    tmp_path: Path,
    entry: Callable[..., int],
    extra: list[str],
    use_sys: bool,
    monkeypatch: MonkeyPatch,
    caplog: LogCaptureFixture,
) -> None:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "jobs:\n  chunk_size: bad\n"
        "resources:\n"
        "  dictionary_dir: dictionary\n"
        "  iuphar_target_csv: dictionary/_IUPHAR/_IUPHAR_target.csv\n"
        "  iuphar_family_csv: dictionary/_IUPHAR/_IUPHAR_family.csv\n"
        "  uniprot_data_dir: uniprot\n"
        "  organism_csv: dictionary/organism.csv\n"
        "  status_csv: dictionary/status.csv\n"
        "  targets_type_csv: dictionary/targets_type.csv\n"
    )
    argv = [*extra, "--config", str(cfg)]
    buf = io.StringIO()
    orig = configure_logger

    def _conf(cfg: LoggerConfig, *a: Any, **k: Any) -> Logger:
        cfg.stream = buf
        return orig(cfg, *a, **k)

    monkeypatch.setattr("library.cli.configure_logger", _conf)
    if use_sys:
        monkeypatch.setattr(sys, "argv", ["prog", *argv])
        rc = entry()
    else:
        rc = entry(argv)
    assert rc != 0

    if caplog.text:
        assert "jobs.chunk_size" in caplog.text
