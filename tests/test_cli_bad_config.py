import sys
import io
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
from library.cli import configure_logger


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
    tmp_path: Path, entry, extra, use_sys, monkeypatch
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

    def _conf(cfg, *a, **k):
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

