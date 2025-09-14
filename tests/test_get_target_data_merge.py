"""Tests for merging IUPHAR classifications with existing target data."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from pytest import MonkeyPatch

import library.target_postprocessing as tp
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
    cfg.target.all.organism_csv = tmp_path / "organism.csv"
    pd.DataFrame({"genus": [], "type": []}).to_csv(
        cfg.target.all.organism_csv, index=False
    )

    merged = gtd.merge_results(combined_df, iuphar_df, cfg)

    assert "ec_number" in merged.columns
    assert "ec_number_x" not in merged.columns
    assert "ec_number_y" not in merged.columns
    assert merged.loc[0, "ec_number"] == "1.1.1.1"
