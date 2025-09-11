from __future__ import annotations

from pathlib import Path

import pandas as pd

import get_assay_data as gas
import get_document_data as gdd
import get_testitem_data as gtdt
import get_target_data as gtd
from library import chembl_library as cl, io
from library.config import Config


def _create_config(tmp_path: Path) -> Path:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "jobs:\n  chunk_size: 10\n"
        "io:\n  csv_sep: ','\n  csv_encoding: utf8\n"
        "log:\n  level: INFO\n"
        "api:\n  timeout_read: 30\n"
        "resources:\n"
        "  dictionary_dir: dictionary\n"
        "  iuphar_target_csv: dictionary/_IUPHAR/_IUPHAR_target.csv\n"
        "  iuphar_family_csv: dictionary/_IUPHAR/_IUPHAR_family.csv\n"
        "  uniprot_data_dir: uniprot\n"
        "  organism_csv: dictionary/organism.csv\n"
    )
    return cfg


def test_assay_timeout_override(tmp_path: Path, monkeypatch) -> None:
    input_csv = tmp_path / "assay.csv"
    input_csv.write_text("assay_chembl_id\na1\n")
    config_path = _create_config(tmp_path)
    called: dict[str, float] = {}
    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["a1"])
    monkeypatch.setattr(
        gas, "apply_config_overrides", lambda a, p, c, mapping=None: Config()
    )

    def fake_get_assays(ids, cfg, chunk_size, timeout):
        called["timeout"] = timeout
        return pd.DataFrame({"assay_chembl_id": ids})

    monkeypatch.setattr(cl, "get_assays", fake_get_assays)
    monkeypatch.setattr(gas.ap, "postprocess_assays", lambda df: df)
    monkeypatch.setattr(io, "write_csv", lambda *a, **k: None)
    monkeypatch.setattr(gas, "analyze_table_quality", lambda df, table_name: None)
    rc = gas.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--timeout",
            "5",
        ]
    )
    assert rc == 0
    assert called["timeout"] == 5


def test_document_timeout_override(tmp_path: Path, monkeypatch) -> None:
    input_csv = tmp_path / "doc.csv"
    input_csv.write_text("chembl_id\nd1\n")
    config_path = _create_config(tmp_path)
    called: dict[str, float] = {}
    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["d1"])
    monkeypatch.setattr(
        gdd, "apply_config_overrides", lambda a, p, c, mapping=None: Config()
    )

    def fake_get_documents(ids, cfg, chunk_size, timeout):
        called["timeout"] = timeout
        return pd.DataFrame({"document_chembl_id": ids})

    monkeypatch.setattr(cl, "get_documents", fake_get_documents, raising=False)
    monkeypatch.setattr(io, "write_csv", lambda *a, **k: None)
    monkeypatch.setattr(gdd, "analyze_table_quality", lambda df, table_name: None)
    rc = gdd.main(
        [
            "chembl",
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--timeout",
            "7",
        ]
    )
    assert rc == 0
    assert called["timeout"] == 7


def test_testitem_timeout_override(tmp_path: Path, monkeypatch) -> None:
    input_csv = tmp_path / "testitem.csv"
    input_csv.write_text("molecule_chembl_id\nt1\n")
    config_path = _create_config(tmp_path)
    called: dict[str, float] = {}
    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["t1"])
    monkeypatch.setattr(
        gtdt, "apply_config_overrides", lambda a, p, c, mapping=None: Config()
    )

    def fake_get_testitem(ids, cfg, chunk_size, timeout):
        called["timeout"] = timeout
        return pd.DataFrame({"molecule_chembl_id": ids})

    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)
    monkeypatch.setattr(gtdt, "add_pubchem_data", lambda df, cfg: df)
    monkeypatch.setattr(io, "write_csv", lambda *a, **k: None)
    monkeypatch.setattr(gtdt, "analyze_table_quality", lambda df, table_name: None)
    rc = gtdt.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--timeout",
            "9",
        ]
    )
    assert rc == 0
    assert called["timeout"] == 9


def test_target_timeout_override(tmp_path: Path, monkeypatch) -> None:
    input_csv = tmp_path / "target.csv"
    input_csv.write_text("chembl_id\nt1\n")
    config_path = _create_config(tmp_path)
    called: dict[str, float] = {}
    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["t1"])
    monkeypatch.setattr(
        gtd, "apply_config_overrides", lambda a, p, c, mapping=None: Config()
    )

    def fake_get_targets(ids, cfg, mapping_cfg, timeout):
        called["timeout"] = timeout
        return pd.DataFrame({"target_chembl_id": ids})

    monkeypatch.setattr(cl, "get_targets", fake_get_targets)
    monkeypatch.setattr(io, "write_csv", lambda *a, **k: None)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda df, table_name: None)
    rc = gtd.main(
        [
            "chembl",
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--timeout",
            "11",
        ]
    )
    assert rc == 0
    assert called["timeout"] == 11
