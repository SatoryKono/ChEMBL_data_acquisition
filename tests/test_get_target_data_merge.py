"""Tests for merging IUPHAR classifications with existing target data."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from pytest import MonkeyPatch

import library.target_postprocessing as tp
from library import protein_classification as pc
from library.config import Config
from scripts import get_target_data as gtd


def test_iuphar_merge_preserves_ec_number(
    tmp_path: Path, cfg: Config, monkeypatch: MonkeyPatch
) -> None:
    """Ensure merging does not duplicate the ``ec_number`` column."""
    combined_df = pd.DataFrame(
        {
            "chembl_id": ["CHEMBL1"],
            "uniprot_id": ["P12345"],
            "ec_number": ["1.1.1.1"],
            "synonyms": ["foo"],
            "gene_name": ["GN1"],
        }
    )
    iuphar_df = pd.DataFrame({"uniprot_id": ["P12345"], "class_a": ["Enzyme"]})

    monkeypatch.setattr(tp, "postprocess_targets", lambda df: df)
    monkeypatch.setattr(tp, "finalise_targets", lambda df, org: df)

    sentinel_classifier = object()
    calls: dict[str, int | bool] = {"from_config": 0, "append": False}

    def fake_classifier_from_config(cfg_obj: Config) -> object:
        calls["from_config"] += 1
        return sentinel_classifier

    def fake_append(df: pd.DataFrame, classifier: object) -> pd.DataFrame:
        assert classifier is sentinel_classifier
        calls["append"] = True
        return df

    monkeypatch.setattr(pc, "classifier_from_config", fake_classifier_from_config)
    monkeypatch.setattr(pc, "append_protein_class_predictions", fake_append)
    cfg.target.all.organism_csv = Path("tests/data/organism_min.csv")

    merged = gtd.merge_results(combined_df, iuphar_df, cfg)

    assert "ec_number" in merged.columns
    assert "ec_number_x" not in merged.columns
    assert "ec_number_y" not in merged.columns
    assert merged.loc[0, "ec_number"] == "1.1.1.1"
    assert calls["from_config"] == 1
    assert calls["append"] is True


def test_merge_results_preserves_existing_type(
    cfg: Config, monkeypatch: MonkeyPatch
) -> None:
    """Existing ``type`` values are retained as ``target_type`` after merge."""

    combined_df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprot_id": ["P1"],
            "genus": ["Homo"],
            "type": ["original"],
        }
    )
    iuphar_df = pd.DataFrame({"uniprot_id": ["P1"], "class": ["enzyme"]})

    def fake_postprocess(df: pd.DataFrame, *, chembl_col: str = "target_chembl_id") -> pd.DataFrame:
        processed = df.rename(columns={"uniprot_id": "uniprotkb_Id"}).copy()
        processed["organism"] = processed["genus"]
        return processed

    sentinel_classifier = object()

    monkeypatch.setattr(tp, "postprocess_targets", fake_postprocess)
    monkeypatch.setattr(pc, "classifier_from_config", lambda cfg_obj: sentinel_classifier)
    monkeypatch.setattr(
        pc,
        "append_protein_class_predictions",
        lambda df, classifier: df,
    )
    cfg.target.all.organism_csv = Path("tests/data/organism_min.csv")

    merged = gtd.merge_results(combined_df, iuphar_df, cfg)

    assert merged.loc[0, "target_type"] == "original"
    assert "type" not in merged.columns
