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

import mapper_main
from library.cli import LoggerConfig, configure_logger
from library.config import Config


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


def test_mapper_main_logs_mapping(tmp_path: Path, monkeypatch: Any) -> None:
    """Mapper CLI emits structured mapping log records."""

    def fake_map(cid: str, cfg: object) -> str | None:
        return "P12345"

    monkeypatch.setattr(mapper_main, "map_chembl_to_uniprot", fake_map)
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
    )
    buffer = io.StringIO()
    configure_logger(LoggerConfig(stream=buffer, level="INFO"))

    base_cfg = Config()
    cfg = SimpleNamespace(
        io=base_cfg.io, uniprot_mapping=base_cfg.uniprot_mapping, to_dict=lambda: {}
    )
    exit_code = mapper_main.run(cfg, args)
    assert exit_code == 0
    records = [json.loads(line) for line in buffer.getvalue().splitlines()]
    mapped = [r for r in records if r["event"] == "mapped"]
    assert mapped
    rec = mapped[0]
    assert rec["chembl_id"] == "CHEMBL1"
    assert rec["uniprot_id"] == "P12345"
