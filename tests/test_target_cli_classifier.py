"""CLI integration tests for target pipeline classifier usage."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from pytest import MonkeyPatch

from library.config import Config
from scripts import get_target_data as gtd


def test_all_parser_accepts_missing_organism_argument(tmp_path: Path) -> None:
    """The ``all`` sub-command should parse without ``--organism-csv``."""

    parser, _log_cfg = gtd.build_parser()
    output = tmp_path / "out.csv"
    args = parser.parse_args(
        [
            "all",
            "--input",
            str(tmp_path / "input.csv"),
            "--output",
            str(output),
        ]
    )
    assert getattr(args, "organism_csv", None) is None


def test_run_all_initialises_classifier(
    monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    """``run_all`` should load and use the IUPHAR classifier."""

    chembl_df = pd.DataFrame({"target_chembl_id": ["CHEMBL1"], "uniprot_id": ["P1"]})
    uniprot_df = pd.DataFrame({"uniprot_id": ["P1"], "original_id": ["P1"]})
    iuphar_df = pd.DataFrame({"uniprot_id": ["P1"], "class": ["enzyme"]})

    monkeypatch.setattr(gtd, "fetch_chembl", lambda *a, **k: chembl_df)
    monkeypatch.setattr(gtd, "fetch_uniprot", lambda *a, **k: uniprot_df)
    monkeypatch.setattr(gtd, "fetch_iuphar", lambda *a, **k: (chembl_df, iuphar_df))
    monkeypatch.setattr(gtd, "validate_and_write", lambda *a, **k: 0)

    sentinel_classifier = object()
    calls: dict[str, int | bool] = {"from_config": 0, "append": False}

    def fake_classifier_from_config(cfg_obj: Config) -> object:
        calls["from_config"] += 1
        return sentinel_classifier

    def fake_append(df: pd.DataFrame, classifier: object) -> pd.DataFrame:
        assert classifier is sentinel_classifier
        calls["append"] = True
        return df

    monkeypatch.setattr(gtd.pc, "classifier_from_config", fake_classifier_from_config)
    monkeypatch.setattr(gtd.pc, "append_protein_class_predictions", fake_append)
    monkeypatch.setattr(gtd.tp, "postprocess_targets", lambda df, **_: df)
    monkeypatch.setattr(gtd.tp, "finalise_targets", lambda df, *a, **k: df)

    args = argparse.Namespace(
        input_csv=tmp_path / "input.csv",
        output_csv=tmp_path / "out.csv",
    )

    rc = gtd.run_all(cfg, args)

    assert rc == 0
    assert calls["from_config"] == 1
    assert calls["append"] is True
