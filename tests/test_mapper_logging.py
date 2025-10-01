import argparse
import importlib
import io
import json
import logging
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd

from library.cli import LoggerConfig, configure_logger
from library.config import Config
from library.utils.cli_tools import mapper_main


def test_mapper_library_has_no_logging_side_effect() -> None:
    root = logging.getLogger()
    for handler in list(root.handlers):
        root.removeHandler(handler)
    assert not root.handlers
    module_name = "library.mapper_library"
    if module_name in sys.modules:
        del sys.modules[module_name]
    importlib.import_module(module_name)
    assert not root.handlers


def test_mapper_main_logs_mapping(
    tmp_path: Path, monkeypatch: Any, cfg: Config
) -> None:
    """Mapper CLI emits structured mapping log records."""

    def fake_map(ids: list[str], cfg: object, *, batch_size: int, rps: float, max_workers: int | None) -> dict[str, str | None]:
        assert ids == ["CHEMBL1"]
        assert batch_size == 1
        assert rps == 1.0
        assert max_workers == 1
        return {"CHEMBL1": "P12345"}

    monkeypatch.setattr(mapper_main, "map_chembl_ids_to_uniprot", fake_map)

    df = pd.DataFrame({"chembl_id": ["CHEMBL1"]})
    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False)
    args = argparse.Namespace(
        input_csv=input_path,
        output_csv=tmp_path / "out.csv",
        column="chembl_id",
        sep=",",
        encoding="utf8",
        key_cols=None,
        chunk_size=1,
        rps=1.0,
        workers=1,
    )
    buffer = io.StringIO()
    configure_logger(LoggerConfig(stream=buffer, level="INFO"))

    cfg_ns = SimpleNamespace(
        io=cfg.io,
        uniprot_mapping=cfg.uniprot_mapping,
        to_dict=lambda: {},
    )
    exit_code = mapper_main.run(cfg_ns, args)
    assert exit_code == 0
    records = [json.loads(line) for line in buffer.getvalue().splitlines()]
    mapped = [r for r in records if r["event"] == "mapped"]
    assert mapped
    rec = mapped[0]
    assert rec["chembl_id"] == "CHEMBL1"
    assert rec["uniprot_id"] == "P12345"


def test_mapper_main_concurrent_preserves_order(
    tmp_path: Path, monkeypatch: Any, cfg: Config
) -> None:
    """Concurrent path retains mapping order and logging semantics."""

    df = pd.DataFrame(
        {"chembl_id": ["CHEMBL1", None, "CHEMBL2", " ", "CHEMBL3"]}
    )
    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False)
    output_path = tmp_path / "out.csv"

    captured: dict[str, Any] = {}

    def fake_map(
        ids: list[str],
        cfg_mapping: object,
        *,
        batch_size: int,
        rps: float,
        max_workers: int | None,
    ) -> dict[str, str | None]:
        captured["ids"] = list(ids)
        captured["batch_size"] = batch_size
        captured["rps"] = rps
        captured["max_workers"] = max_workers
        return {
            "CHEMBL2": "P222",
            "CHEMBL1": "P111",
            "CHEMBL3": None,
        }

    monkeypatch.setattr(mapper_main, "map_chembl_ids_to_uniprot", fake_map)

    written: dict[str, Any] = {}

    def fake_write_csv(
        df_out: pd.DataFrame,
        path: Path,
        *,
        cfg: Config,
        sep: str,
        encoding: str,
        key_cols: Any,
    ) -> Path:
        written["df"] = df_out.copy()
        written["path"] = path
        written["key_cols"] = key_cols
        written["chembl_ids"] = df_out["chembl_id"].tolist()
        return path

    monkeypatch.setattr(mapper_main.io, "write_csv", fake_write_csv)

    args = argparse.Namespace(
        input_csv=input_path,
        output_csv=output_path,
        column="chembl_id",
        sep=",",
        encoding="utf8",
        key_cols=None,
        chunk_size=2,
        rps=5.0,
        workers=3,
    )

    buffer = io.StringIO()
    configure_logger(LoggerConfig(stream=buffer, level="INFO"))

    cfg_ns = SimpleNamespace(
        io=cfg.io,
        uniprot_mapping=cfg.uniprot_mapping,
        to_dict=lambda: {},
    )

    exit_code = mapper_main.run(cfg_ns, args)
    assert exit_code == 0

    assert captured["ids"] == ["CHEMBL1", "CHEMBL2", "CHEMBL3"]
    assert captured["batch_size"] == 2
    assert captured["rps"] == 5.0
    assert captured["max_workers"] == 3

    output_df = written["df"]
    chembl_ids = [
        value if isinstance(value, str) and value.strip() else None
        for value in written["chembl_ids"]
    ]
    assert chembl_ids == ["CHEMBL1", None, "CHEMBL2", "CHEMBL3"]
    mapped_values = [
        value if isinstance(value, str) and value else None
        for value in output_df["mapping_uniprot_id"].tolist()
    ]
    assert mapped_values == ["P111", None, "P222", None]
    assert written["path"] == output_path
    assert written["key_cols"] is None

    records = [json.loads(line) for line in buffer.getvalue().splitlines() if line]
    status_events = [
        (record["event"], record.get("chembl_id"))
        for record in records
        if record["event"] in {"mapped", "uniprot_id_missing"}
    ]
    assert status_events == [
        ("mapped", "CHEMBL1"),
        ("mapped", "CHEMBL2"),
        ("uniprot_id_missing", "CHEMBL3"),
    ]
