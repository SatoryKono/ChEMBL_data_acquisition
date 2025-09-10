import logging
import sys
from pathlib import Path

import get_activity_data as gad
import get_assay_data as gas
import get_document_data as gdd
import get_testitem_data as gtdt
import mapper_main as mapper
import table_quality_main as tqm
import get_input_initialisation as gii
import get_target_data as gtd
import get_document_type as gdoctype
import pytest


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
    tmp_path: Path, caplog: pytest.LogCaptureFixture, entry, extra, use_sys, monkeypatch
) -> None:
    cfg = tmp_path / "config.yaml"
    cfg.write_text("jobs:\n  chunk_size: bad\n")
    argv = ["--config", str(cfg), *extra]
    with caplog.at_level(logging.ERROR):
        if use_sys:
            monkeypatch.setattr(sys, "argv", ["prog", *argv])
            rc = entry()
        else:
            rc = entry(argv)
    assert rc != 0
    assert "jobs.chunk_size" in caplog.text
