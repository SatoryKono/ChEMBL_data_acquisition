"""Tests for helper utilities in :mod:`scripts.get_target_data`."""

import argparse
from pathlib import Path

import pandas as pd
import pytest

import scripts.get_target_data as gtd
from library.config import Config


def test_pipe_merge_deduplicates_and_sorts() -> None:
    values = ["b|a", "C|a", None, "", "b"]
    assert gtd._pipe_merge(values) == "C|a|b"


def test_pipe_merge_handles_empty() -> None:
    assert gtd._pipe_merge([None, "", " "]) == ""


def test_first_token_extracts_head() -> None:
    assert gtd._first_token("a|b|c") == "a"
    assert gtd._first_token(None) == ""
    assert gtd._first_token("") == ""


def test_save_snapshot_creates_unique_files(tmp_path: Path, cfg: Config) -> None:
    df = pd.DataFrame({"a": [1]})
    base = tmp_path / "out.csv"
    gtd._save_snapshot(df, base, "step", cfg)
    gtd._save_snapshot(df, base, "step", cfg)
    files = sorted(tmp_path.glob("out_step_*.csv"))
    assert len(files) == 2


def test_run_all_saves_snapshots(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``run_all`` should snapshot data after each stage."""

    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding=cfg.io.csv_encoding)
    output_csv = tmp_path / "final.csv"

    def fake_run_chembl(cfg: Config, args: argparse.Namespace) -> int:
        df = pd.DataFrame(
            {
                "target_chembl_id": ["CHEMBL1"],
                "uniprot_id": ["P12345"],
                "gene": ["GN1"],
                "component_description": ["desc"],
                "pref_name": ["name"],
                "chembl_alternative_name": ["alt"],
                "ec_numbers": ["1.1.1.1"],
                "reaction_ec_numbers": ["2.2.2.2"],
            }
        )
        df.to_csv(
            args.output_csv,
            index=False,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
        return 0

    def fake_run_uniprot(cfg: Config, args: argparse.Namespace) -> int:
        df = pd.DataFrame(
            {
                "uniprot_id": ["P12345"],
                "names": ["prot"],
                "secondaryAccessionNames": ["sec"],
                "ec_numbers": ["1.1.1.1"],
                "reaction_ec_numbers": ["2.2.2.2"],
            }
        )
        df.to_csv(
            args.output_csv,
            index=False,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
        return 0

    def fake_run_iuphar(cfg: Config, args: argparse.Namespace) -> int:
        df = pd.DataFrame({"uniprot_id": ["P12345"], "class_a": ["Enzyme"]})
        df.to_csv(
            args.output_csv,
            index=False,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
        return 0

    steps: list[str] = []

    def fake_snapshot(
        df: pd.DataFrame, base: Path, step: str, cfg: Config
    ) -> Path:  # pragma: no cover - simple recorder
        steps.append(step)
        return base

    monkeypatch.setattr(gtd, "run_chembl", fake_run_chembl)
    monkeypatch.setattr(gtd, "run_uniprot", fake_run_uniprot)
    monkeypatch.setattr(gtd, "run_iuphar", fake_run_iuphar)
    monkeypatch.setattr(gtd, "_save_snapshot", fake_snapshot)
    monkeypatch.setattr(gtd.tp, "postprocess_targets", lambda df: df)
    monkeypatch.setattr(gtd.tp, "finalise_targets", lambda df, _: df)

    args = argparse.Namespace(input_csv=input_csv, output_csv=output_csv)
    rc = gtd.run_all(cfg, args)
    assert rc == 0
    assert steps == ["chembl", "uniprot", "iuphar"]
