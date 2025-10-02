"""Tests for validate_and_write utility."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from pytest import MonkeyPatch

from library.config import Config
from scripts import get_target_data as gtd


def test_validate_and_write(
    tmp_path: Path, cfg: Config, monkeypatch: MonkeyPatch
) -> None:
    df = pd.DataFrame({"target_chembl_id": ["C1"]})
    out = tmp_path / "out.csv"
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *a, **k: None)
    exit_code = gtd.validate_and_write(df, out, cfg)
    assert exit_code == 0
    assert out.exists()


def test_validate_and_write_cleans_identifier_placeholders(
    tmp_path: Path, cfg: Config, monkeypatch: MonkeyPatch
) -> None:
    df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprot_id_primary": [" - "],
            "gene_symbol": [None],
        }
    )
    output = tmp_path / "targets.csv"

    def _write_csv(
        frame: pd.DataFrame,
        path: Path,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
        col_order: list[str] | None = None,
        key_cols: list[str] | None = None,
        chunksize: int | None = None,
    ) -> Path:
        frame.to_csv(path, index=False, sep=sep or cfg.io.csv_sep, encoding=encoding or cfg.io.csv_encoding)
        return Path(path)

    monkeypatch.setattr(gtd, "normalize_targets", lambda frame: frame)
    monkeypatch.setattr(gtd, "add_pipeline_metadata", lambda frame: frame)
    monkeypatch.setattr(
        gtd,
        "_prepare_targets_for_schema",
        lambda frame: (frame, set(), set()),
    )
    monkeypatch.setattr(
        gtd.TargetsSchema,
        "validate",
        staticmethod(lambda frame, lazy=True: frame),
    )
    monkeypatch.setattr(gtd.io, "write_csv", _write_csv)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *args, **kwargs: None)

    exit_code = gtd.validate_and_write(df, output, cfg)

    assert exit_code == 0
    result = pd.read_csv(output, dtype=str, keep_default_na=False)
    assert result.loc[0, "uniprot_id_primary"] == ""
    assert result.loc[0, "gene_symbol"] == "-"
