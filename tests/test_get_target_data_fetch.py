"""Tests for data fetching utilities in get_target_data."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from pytest import MonkeyPatch

from library.config import Config
from scripts import get_target_data as gtd


def test_fetch_chembl(monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config) -> None:
    out = tmp_path / "chembl.csv"
    inp = tmp_path / "in.csv"

    def fake_run_chembl(cfg: Config, args: argparse.Namespace) -> int:
        pd.DataFrame({"target_chembl_id": ["C1"], "uniprot_id": ["P1"]}).to_csv(
            args.output_csv, index=False
        )
        return 0

    monkeypatch.setattr(gtd, "run_chembl", fake_run_chembl)
    df = gtd.fetch_chembl(cfg, inp, out)
    assert df.loc[0, "uniprot_id"] == "P1"


def test_fetch_uniprot(monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config) -> None:
    chembl_df = pd.DataFrame({"uniprot_id": ["P12345"]})
    out = tmp_path / "uniprot.csv"

    def fake_run_uniprot(cfg: Config, args: argparse.Namespace) -> int:
        pd.DataFrame({"uniprot_id": ["P12345"], "names": ["Foo"]}).to_csv(
            args.output_csv, index=False
        )
        return 0

    monkeypatch.setattr(gtd, "run_uniprot", fake_run_uniprot)
    df = gtd.fetch_uniprot(cfg, chembl_df, out)
    assert list(df["original_id"]) == ["P12345"]


def test_run_uniprot_initialises_session(
    monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("uniprot_id\nP12345\n", encoding=cfg.io.csv_encoding)
    output_csv = tmp_path / "uniprot.csv"
    cfg.target.uniprot.data_dir = tmp_path

    called: dict[str, object] = {}

    def fake_init_session(api: object, retry: object) -> None:
        called["init"] = (api, retry)

    def fake_process(
        input_csv: str,
        output_csv: str,
        data_dir: Path | str | None = None,
        *,
        cfg: object,
        sep: str = ",",
        encoding: str = "utf-8",
    ) -> None:
        called["cfg"] = cfg
        pd.DataFrame({"uniprot_id": ["P12345"], "names": ["Foo"]}).to_csv(
            output_csv, index=False
        )

    monkeypatch.setattr(gtd.uu, "init_session", fake_init_session)
    monkeypatch.setattr(gtd.uu, "process", fake_process)

    def fake_write_csv(
        df: pd.DataFrame,
        path: Path,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
        **__: object,
    ) -> Path:
        return path

    monkeypatch.setattr(gtd.io, "write_csv", fake_write_csv)

    args = argparse.Namespace(input_csv=input_csv, output_csv=output_csv)
    rc = gtd.run_uniprot(cfg, args)

    assert rc == 0
    assert called["init"] == (cfg.api, cfg.retry)
    assert called["cfg"] is cfg.uniprot


def test_fetch_iuphar(monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config) -> None:
    chembl_df = pd.DataFrame(
        {
            "target_chembl_id": ["C1"],
            "uniprot_id": ["P1"],
            "pref_name": ["pref"],
            "component_description": ["desc"],
            "gene": ["GENE"],
            "chembl_alternative_name": ["alt"],
            "ec_numbers": ["1.1.1.1"],
            "reaction_ec_numbers": ["2.2.2.2"],
        }
    )
    uniprot_df = pd.DataFrame(
        {
            "uniprot_id": ["P1"],
            "original_id": ["P1"],
            "names": ["name"],
            "secondaryAccessionNames": ["sec"],
            "ec_numbers": ["3.3.3.3"],
            "reaction_ec_numbers": ["4.4.4.4"],
        }
    )
    out = tmp_path / "iuphar.csv"

    def fake_run_iuphar(cfg: Config, args: argparse.Namespace) -> int:
        pd.DataFrame({"uniprot_id": ["P1"], "IUPHAR_class": ["Enzyme"]}).to_csv(
            args.output_csv, index=False
        )
        return 0

    monkeypatch.setattr(gtd, "run_iuphar", fake_run_iuphar)
    combined_df, iuphar_df = gtd.fetch_iuphar(cfg, chembl_df, uniprot_df, out)
    assert "synonyms" in combined_df.columns
    assert iuphar_df.loc[0, "IUPHAR_class"] == "Enzyme"
