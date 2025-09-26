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
from library.utils.cli_tools import get_document_type as gdoctype
from library.utils.cli_tools import get_input_initialisation as gii
from scripts import get_target_data as gtd
from scripts import get_testitem_data as gtdt
from library.utils.cli_tools import mapper_main as mapper
from library.utils.cli_tools import table_quality_main as tqm

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
        "sources:\n"
        "  chembl:\n"
        "    pipelines:\n"
        "      activity:\n"
        "        batch_size: bad\n"
        "local:\n"
        "  resources:\n"
        "    dictionary_dir: dictionary\n"
        "    iuphar_target_csv: dictionary/_IUPHAR/_IUPHAR_target.csv\n"
        "    iuphar_family_csv: dictionary/_IUPHAR/_IUPHAR_family.csv\n"
        "    uniprot_data_dir: uniprot\n"
        "    targets_type_csv: dictionary/targets_type.csv\n"
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
        assert "sources.chembl.pipelines.activity.batch_size" in caplog.text


@pytest.mark.parametrize("entry, extra, use_sys", CLIS)
def test_unknown_key_config_exits(
    tmp_path: Path,
    entry: Callable[..., int],
    extra: list[str],
    use_sys: bool,
    monkeypatch: MonkeyPatch,
) -> None:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "sources:\n"
        "  chembl:\n"
        "    pipelines:\n"
        "      activity:\n"
        "        batch_size: 5\n"
        "local:\n"
        "  resources:\n"
        "    dictionary_dir: dictionary\n"
        "    iuphar_target_csv: dictionary/_IUPHAR/_IUPHAR_target.csv\n"
        "    iuphar_family_csv: dictionary/_IUPHAR/_IUPHAR_family.csv\n"
        "    uniprot_data_dir: uniprot\n"
        "    targets_type_csv: dictionary/targets_type.csv\n"
        "unknown: 1\n"
    )
    argv = [*extra, "--config", str(cfg)]
    buf = io.StringIO()
    orig = configure_logger

    def _conf(cfg: LoggerConfig, *a: Any, **k: Any) -> Logger:
        cfg.stream = buf
        return orig(cfg, *a, **k)

    monkeypatch.setattr("library.cli.configure_logger", _conf)
    module = sys.modules.get(entry.__module__)
    if module and hasattr(module, "configure_logger"):
        monkeypatch.setattr(f"{entry.__module__}.configure_logger", _conf)

    if use_sys:
        monkeypatch.setattr(sys, "argv", ["prog", *argv])
        rc = entry()
    else:
        rc = entry(argv)
    assert rc != 0
    assert "Unknown configuration key" in buf.getvalue()


def test_negative_limit_in_config_exits(
    tmp_path: Path,
    monkeypatch: MonkeyPatch,
) -> None:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "sources:\n"
        "  chembl:\n"
        "    pipelines:\n"
        "      activity:\n"
        "        batch_size: 5\n"
        "        limit: -1\n"
        "local:\n"
        "  resources:\n"
        "    dictionary_dir: dictionary\n"
        "    iuphar_target_csv: dictionary/_IUPHAR/_IUPHAR_target.csv\n"
        "    iuphar_family_csv: dictionary/_IUPHAR/_IUPHAR_family.csv\n"
        "    uniprot_data_dir: uniprot\n"
        "    targets_type_csv: dictionary/targets_type.csv\n"
    )
    buf = io.StringIO()
    orig = configure_logger

    def _conf(cfg: LoggerConfig, *a: Any, **k: Any) -> Logger:
        cfg.stream = buf
        return orig(cfg, *a, **k)

    monkeypatch.setattr("library.cli.configure_logger", _conf)
    monkeypatch.setattr("scripts.get_activity_data.configure_logger", _conf)
    rc = gad.main(["--config", str(cfg)])
    assert rc != 0
    assert "sources.chembl.pipelines.activity.limit" in buf.getvalue()


def test_doc_type_negative_limit_in_config_exits(
    tmp_path: Path,
    monkeypatch: MonkeyPatch,
) -> None:
    cfg = tmp_path / "config.yaml"
    cfg.write_text("system:\n  doc_type:\n    limit: -1\n")
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("chembl_id\nCHEMBL1\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"

    buf = io.StringIO()
    orig = configure_logger

    def _conf(cfg: LoggerConfig, *a: Any, **k: Any) -> Logger:
        cfg.stream = buf
        return orig(cfg, *a, **k)

    monkeypatch.setattr("library.cli.configure_logger", _conf)
    monkeypatch.setattr("library.utils.cli_tools.get_document_type.configure_logger", _conf)

    rc = gdoctype.main(
        [
            "--config",
            str(cfg),
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
        ]
    )

    assert rc != 0
    assert "system.doc_type.limit" in buf.getvalue()
