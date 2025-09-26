from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pytest

from library import chembl_library as cl
from library import io
from library.config import Config
from schemas import TestitemsSchema
from scripts import get_testitem_data as gtd


def test_run_chembl_column_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Ensure schema columns precede alphabetically sorted extras."""
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1"]))

    df = pd.DataFrame(
        [
            {
                "molecule_chembl_id": "CHEMBL1",
                "parent_molecule_id": "CHEMBL0",
                "molecule_type": "Small molecule",
                "salt_chembl_id": "CHEMBL1-SALT",
                "chirality": 1,
                "extra_b": 2,
                "extra_a": 1,
            }
        ]
    )

    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: df)
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda df, cfg: df)
    monkeypatch.setattr(
        gtd,
        "attach_parent_molecule_ids",
        lambda frame, **kwargs: (
            frame,
            gtd.ParentLookupStats(
                source=gtd.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
            ),
        ),
    )
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda p: "deadbeef")

    captured: dict[str, list[str]] = {}

    def fake_write_csv(
        df: pd.DataFrame,
        output: Path,
        *,
        cfg: Config,
        key_cols: list[str] | None = None,
        col_order: list[str] | None = None,
        **__: object,
    ) -> Path:
        captured["col_order"] = list(col_order or [])
        return output

    monkeypatch.setattr(io, "write_csv", fake_write_csv)

    rc = gtd.run_chembl(cfg, args)
    assert rc == 0

    available = set(df.columns) | {"pipeline_version", "timestamp_utc"}
    expected_head = [c for c in TestitemsSchema.columns if c in available]
    expected_tail = sorted(available - set(TestitemsSchema.columns))
    assert captured["col_order"] == expected_head + expected_tail


def test_run_chembl_initialises_pubchem_session(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n", encoding=cfg.io.csv_encoding)

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1"]))

    df = pd.DataFrame(
        {"molecule_chembl_id": ["CHEMBL1"], "molecule_type": ["Small molecule"]}
    )
    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: df)
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, pubchem_cfg: frame)
    monkeypatch.setattr(
        gtd,
        "attach_parent_molecule_ids",
        lambda frame, **kwargs: (
            frame,
            gtd.ParentLookupStats(
                source=gtd.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
            ),
        ),
    )

    captured: dict[str, object] = {}

    def fake_init_session(api: object, retry: object) -> None:
        captured["init"] = (api, retry)

    monkeypatch.setattr(gtd.pl, "init_session", fake_init_session)
    monkeypatch.setattr(
        io,
        "write_csv",
        lambda df, path, *, cfg, key_cols=None, col_order=None, **__: path,
    )
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "deadbeef")

    rc = gtd.run_chembl(cfg, args)

    assert rc == 0
    assert captured["init"] == (cfg.api, cfg.retry)


def test_ensure_no_parant_column_detects_typo() -> None:
    """Utility rejects frames containing the typo column."""

    df = pd.DataFrame({"parant_molecule_id": ["CHEMBL0"]})

    with pytest.raises(ValueError):
        gtd.ensure_no_parant_column(df)


def test_run_chembl_fails_on_parant_column(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """``run_chembl`` aborts early when the typo column is present."""

    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1"]))

    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "parant_molecule_id": ["CHEMBL0"],
        }
    )

    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: df)
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, pubchem_cfg: frame)

    def fail_write(*_: object, **__: object) -> Path:
        raise AssertionError("write_csv should not be called when validation fails")

    monkeypatch.setattr(io, "write_csv", fail_write)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "deadbeef")

    rc = gtd.run_chembl(cfg, args)

    assert rc == 1
