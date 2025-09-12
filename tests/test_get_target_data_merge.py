"""Tests for merging IUPHAR classifications with existing target data."""

from __future__ import annotations

import pandas as pd


def test_iuphar_merge_preserves_ec_number() -> None:
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
    iuphar_df = pd.DataFrame(
        {
            "uniprot_id": ["P12345"],
            "class_a": ["Enzyme"],
            "synonyms": ["bar"],
            "ec_number": ["2.2.2.2"],
            "gene_name": ["GN1"],
        }
    )

    existing_cols = set(combined_df.columns)
    classification_cols = [c for c in iuphar_df.columns if c not in existing_cols]
    iuphar_df = iuphar_df[["uniprot_id", *classification_cols]].copy()

    merged = combined_df.merge(iuphar_df, on="uniprot_id", how="left")

    assert "ec_number" in merged.columns
    assert "ec_number_x" not in merged.columns
    assert "ec_number_y" not in merged.columns
    assert merged.loc[0, "ec_number"] == "1.1.1.1"
