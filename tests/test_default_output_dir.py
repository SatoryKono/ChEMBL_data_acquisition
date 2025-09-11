from __future__ import annotations

from pathlib import Path
import argparse

import pandas as pd
import pytest

from library import io
from library.config import Config, IoCfg
import mapper_main


def test_default_output_path_uses_output_dir(tmp_path: Path) -> None:
    cfg = IoCfg(output_dir=tmp_path)
    result = io.default_output_path(tmp_path / "input.csv", cfg)
    assert result.parent == tmp_path
    assert result.name.startswith("output_input_")


def test_mapper_run_defaults_to_io_output_dir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_path = tmp_path / "data.csv"
    input_path.write_text("chembl_id\nCHEMBL1\n")

    cfg = Config()
    cfg.io.output_dir = tmp_path / "results"

    args = argparse.Namespace(
        input_csv=input_path,
        output_csv=None,
        column="chembl_id",
        sep=",",
        encoding="utf8",
    )

    monkeypatch.setattr(mapper_main, "map_chembl_to_uniprot", lambda _i, _cfg: "P12345")
    monkeypatch.setattr(io, "_write_meta", lambda *_a, **_k: None)

    rc = mapper_main.run(cfg, args)
    assert rc == 0

    expected = io.default_output_path(input_path, cfg.io)
    assert expected.exists()
    df = pd.read_csv(expected)
    assert "mapping_uniprot_id" in df.columns
