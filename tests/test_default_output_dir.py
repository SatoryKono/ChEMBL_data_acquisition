from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pytest

from library import io
from library.config import Config, IoCfg
from library.utils.cli_tools import mapper_main


def test_default_output_path_uses_output_dir(tmp_path: Path) -> None:
    cfg = IoCfg(output_dir=tmp_path)
    result = io.default_output_path(tmp_path / "input.csv", cfg)
    assert result.parent == tmp_path
    assert result.name.startswith("output.input_")


def test_mapper_run_defaults_to_io_output_dir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_path = tmp_path / "data.csv"
    input_path.write_text("chembl_id\nCHEMBL1\n")

    cfg.io.output_dir = tmp_path / "results"

    args = argparse.Namespace(
        input_csv=input_path,
        output_csv=None,
        column="chembl_id",
        sep=",",
        encoding="utf8",
        key_cols=None,
        chunk_size=1,
        rps=1.0,
        workers=1,
    )

    monkeypatch.setattr(
        mapper_main,
        "map_chembl_ids_to_uniprot",
        lambda ids, _cfg, *, batch_size, rps, max_workers: {ids[0]: "P12345"},
    )
    monkeypatch.setattr(io, "write_meta_yaml", lambda *_a, **_k: None)

    rc = mapper_main.run(cfg, args)
    assert rc == 0

    expected = io.default_output_path(input_path, cfg.io)
    assert expected.exists()
    df = pd.read_csv(expected)
    assert "mapping_uniprot_id" in df.columns
