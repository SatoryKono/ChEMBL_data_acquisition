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
    monkeypatch.setattr(
        gtd, "load_parent_catalog", lambda **__: {"CHEMBL1": "CHEMBL1_PARENT"}
    )
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



def test_run_chembl_merges_parent_catalog(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\nCHEMBL2\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1", "CHEMBL2"]))

    source = pd.DataFrame(
        [
            {"molecule_chembl_id": "CHEMBL1", "parent_molecule_chembl_id": None},
            {
                "molecule_chembl_id": "CHEMBL2",
                "parent_molecule_chembl_id": "CHEMBL2_EXISTING",
            },
        ]
    )

    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: source.copy())
    monkeypatch.setattr(
        gtd,
        "load_parent_catalog",
        lambda **__: {"CHEMBL1": "CHEMBL1_PARENT", "CHEMBL2": "CHEMBL2_PARENT"},
    )
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, _: frame)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "deadbeef")

    captured_df: pd.DataFrame | None = None

    def fake_write_csv(
        df: pd.DataFrame,
        output: Path,
        *,
        cfg: Config,
        key_cols: list[str] | None = None,
        col_order: list[str] | None = None,
        **__: object,
    ) -> Path:
        nonlocal captured_df
        captured_df = df.copy()
        return output

    monkeypatch.setattr(io, "write_csv", fake_write_csv)

    rc = gtd.run_chembl(cfg, args)
    assert rc == 0
    assert captured_df is not None
    assert captured_df["parent_molecule_chembl_id"].tolist() == [
        "CHEMBL1_PARENT",
        "CHEMBL2_EXISTING",
    ]


def test_run_chembl_parent_catalog_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:

    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1"]))

    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: pd.DataFrame())
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, _: frame)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "deadbeef")

    monkeypatch.setattr(
        gtd,
        "load_parent_catalog",
        lambda **__: (_ for _ in ()).throw(
            ValueError("missing columns: parant_molecule_id")
        ),
    )

    called = False

    def fake_write_csv(*args: object, **kwargs: object) -> Path:  # pragma: no cover
        nonlocal called
        called = True
        return Path("unused.csv")

    monkeypatch.setattr(io, "write_csv", fake_write_csv)

    rc = gtd.run_chembl(cfg, args)
    assert rc == 1
    assert not called

