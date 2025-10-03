from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pytest

from library import io
from library.config import Config, IoCfg
from library.io import paths as io_paths
from library.utils.cli_tools import mapper_main


UTC_FREEZE_TIME = datetime(2024, 1, 2, 3, 4, 5, tzinfo=timezone.utc)


def _freeze_datetime(monkeypatch: pytest.MonkeyPatch, *, frozen: datetime) -> None:
    class FrozenDateTime:
        @staticmethod
        def now(tz: timezone | None = None) -> datetime:
            if tz is None:
                return frozen
            return frozen.astimezone(tz)

    monkeypatch.setattr(io_paths, "datetime", FrozenDateTime)


def test_default_output_path_uses_output_dir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _freeze_datetime(monkeypatch, frozen=UTC_FREEZE_TIME)
    cfg = IoCfg(output_dir=tmp_path)
    result = io.default_output_path(tmp_path / "input.csv", cfg)
    assert result.parent == tmp_path
    assert result.name == f"output.input_{UTC_FREEZE_TIME.strftime('%Y%m%d')}.csv"


def test_mapper_run_defaults_to_io_output_dir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    _freeze_datetime(monkeypatch, frozen=UTC_FREEZE_TIME)
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
    monkeypatch.setattr(mapper_main, "mapping_failed", False, raising=False)

    rc = mapper_main.run(cfg, args)
    assert rc == 0

    expected = io.default_output_path(input_path, cfg.io)
    assert expected.exists()
    assert expected.name == f"output.data_{UTC_FREEZE_TIME.strftime('%Y%m%d')}.csv"
    df = pd.read_csv(expected)
    assert "mapping_uniprot_id" in df.columns
