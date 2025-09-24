from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd
from pytest import MonkeyPatch

from library import chembl_library as cl
from library import io
from library.config import Config
from scripts import get_assay_data as gas
from scripts import get_document_data as gdd
from scripts import get_target_data as gtd
from scripts import get_testitem_data as gtdt


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
        "  targets_type_csv: dictionary/targets_type.csv\n"
    )
    return cfg


def test_assay_timeout_override(
    tmp_path: Path, monkeypatch: MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "assay.csv"
    input_csv.write_text("assay_chembl_id\na1\n")
    config_path = _create_config(tmp_path)
    called: dict[str, float] = {}
    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["a1"])

    def fake_apply(
        a: Any, p: Any, c: Config, mapping: dict[str, str] | None = None
    ) -> Config:
        cfg.assay.timeout = float(a.timeout)
        return cfg

    monkeypatch.setattr(gas, "apply_config_overrides", fake_apply)

    def fake_get_assays(
        ids: Sequence[str],
        cfg: Config,
        client: Any,
        chunk_size: int,
        timeout: float,
    ) -> pd.DataFrame:
        data = list(ids)
        called["timeout"] = timeout
        return pd.DataFrame({"assay_chembl_id": data})

    monkeypatch.setattr(cl, "get_assays", fake_get_assays)
    monkeypatch.setattr(gas.ap, "postprocess_assays", lambda df: df)
    monkeypatch.setattr(
        io,
        "write_csv",
        lambda df, path, *, cfg, sep=None, encoding=None, **k: path,
    )
    monkeypatch.setattr(gas, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gas, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(gas, "write_meta_yaml", lambda **__: None)
    rc = gas.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--output",
            str(tmp_path / "out.csv"),
            "--timeout",
            "5",
        ]
    )
    assert rc == 0
    assert called["timeout"] == 5


def test_document_timeout_override(
    tmp_path: Path, monkeypatch: MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "doc.csv"
    input_csv.write_text("document_chembl_id\nd1\n")
    config_path = _create_config(tmp_path)
    called: dict[str, float] = {}

    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["d1"])

    def fake_apply_doc(
        a: Any,
        p: Any,
        c: Config,
        mapping: Any | None = None,
        **kwargs: Any,
    ) -> Config:
        cfg.document.chembl.timeout = float(a.timeout)
        cfg.api.user_agent = "test/0.1 (mailto:test@example.com)"
        return cfg

    monkeypatch.setattr(gdd, "apply_config_overrides", fake_apply_doc)

    def fake_get_documents(
        ids: Sequence[str],
        cfg: Config,
        client: Any,
        chunk_size: int,
        timeout: float,
    ) -> pd.DataFrame:
        data = list(ids)
        called["timeout"] = timeout
        return pd.DataFrame({"document_chembl_id": data})

    monkeypatch.setattr(cl, "get_documents", fake_get_documents, raising=False)
    monkeypatch.setattr(
        io,
        "write_csv",
        lambda df, path, *, cfg, sep=None, encoding=None, **k: path,
    )
    monkeypatch.setattr(gdd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gdd, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)
    rc = gdd.main(
        [
            "chembl",
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--output",
            str(tmp_path / "out.csv"),
            "--timeout",
            "7",
        ]
    )
    assert rc == 0
    assert called["timeout"] == 7
    assert cfg.document.chembl.timeout == 7
    assert cfg.api.timeout_read == 7


def test_testitem_timeout_override(
    tmp_path: Path, monkeypatch: MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "testitem.csv"
    input_csv.write_text("molecule_chembl_id\nt1\n")
    config_path = _create_config(tmp_path)
    called: dict[str, float] = {}

    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["t1"])

    def fake_apply_test(
        a: Any, p: Any, c: Config, mapping: Any | None = None
    ) -> Config:
        cfg.testitem.timeout = float(a.timeout)
        return cfg

    monkeypatch.setattr(gtdt, "apply_config_overrides", fake_apply_test)

    def fake_get_testitem(
        ids: Sequence[str],
        cfg: Config,
        client: Any,
        chunk_size: int,
        timeout: float,
    ) -> pd.DataFrame:
        data = list(ids)
        called["timeout"] = timeout
        return pd.DataFrame({"molecule_chembl_id": data})

    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)
    monkeypatch.setattr(gtdt, "add_pubchem_data", lambda df, cfg: df)
    monkeypatch.setattr(
        io,
        "write_csv",
        lambda df, path, *, cfg, sep=None, encoding=None, **k: path,
    )
    monkeypatch.setattr(gtdt, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gtdt, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(gtdt, "write_meta_yaml", lambda **__: None)
    rc = gtdt.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--output",
            str(tmp_path / "out.csv"),
            "--timeout",
            "9",
        ]
    )
    assert rc == 0
    assert called["timeout"] == 9


def test_target_timeout_override(
    tmp_path: Path, monkeypatch: MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "target.csv"
    input_csv.write_text("target_chembl_id\nt1\n")
    config_path = _create_config(tmp_path)
    called: dict[str, float] = {}

    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["t1"])

    def fake_apply_target(
        a: Any,
        p: Any,
        c: Config,
        mapping: Any | None = None,
        **kwargs: Any,
    ) -> Config:
        cfg.target.chembl.timeout = float(a.timeout)
        return cfg

    monkeypatch.setattr(gtd, "apply_config_overrides", fake_apply_target)

    def fake_get_targets(
        ids: Sequence[str],
        cfg: Config,
        client: Any,
        mapping_cfg: Any,
        chunk_size: int,
        timeout: float,
    ) -> pd.DataFrame:
        data = list(ids)
        called["timeout"] = timeout
        return pd.DataFrame({"target_chembl_id": data})

    monkeypatch.setattr(cl, "get_targets", fake_get_targets)
    monkeypatch.setattr(
        io,
        "write_csv",
        lambda df, path, *, cfg, sep=None, encoding=None, **k: path,
    )
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **__: None)
    rc = gtd.main(
        [
            "chembl",
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--output",
            str(tmp_path / "out.csv"),
            "--timeout",
            "11",
        ]
    )
    assert rc == 0
    assert called["timeout"] == 11


def test_target_chunk_size_override(
    tmp_path: Path, monkeypatch: MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "target.csv"
    input_csv.write_text("target_chembl_id\nt1\n")
    config_path = _create_config(tmp_path)
    called: dict[str, int] = {}

    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["t1"])

    def fake_apply_chunk(
        a: Any,
        p: Any,
        c: Config,
        mapping: Any | None = None,
        **kwargs: Any,
    ) -> Config:
        cfg.target.chembl.chunk_size = int(a.chunk_size)
        return cfg

    monkeypatch.setattr(gtd, "apply_config_overrides", fake_apply_chunk)

    def fake_get_targets(
        ids: Sequence[str],
        cfg: Config,
        client: Any,
        mapping_cfg: Any,
        chunk_size: int,
        timeout: float,
    ) -> pd.DataFrame:
        data = list(ids)
        called["chunk_size"] = chunk_size
        return pd.DataFrame({"target_chembl_id": data})

    monkeypatch.setattr(cl, "get_targets", fake_get_targets)
    monkeypatch.setattr(
        io,
        "write_csv",
        lambda df, path, *, cfg, sep=None, encoding=None, **k: path,
    )
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **__: None)
    rc = gtd.main(
        [
            "chembl",
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--output",
            str(tmp_path / "out.csv"),
            "--chunk-size",
            "7",
        ]
    )
    assert rc == 0
    assert called["chunk_size"] == 7
