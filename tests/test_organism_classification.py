"""Tests for :mod:`library.organism_classification`."""

from __future__ import annotations

import pandas as pd

from library import organism_classification as oc


def test_add_cellularity_smart_handles_missing_columns() -> None:
    df = pd.DataFrame({"something": ["value"]})

    result = oc.add_cellularity_smart(df)

    assert result.empty


def test_add_cellularity_smart_detects_viruses() -> None:
    df = pd.DataFrame({
        "genus": ["Lentivirus", "Sample"],
        "superkingdom": ["Viruses", ""],
    })

    result = oc.add_cellularity_smart(df)
    row = result[result["genus"] == "Lentivirus"].iloc[0]
    assert row["type"] == "Viruses"


def test_add_cellularity_smart_detects_bacteria() -> None:
    df = pd.DataFrame(
        {
            "genus": ["Escherichia"],
            "superkingdom": ["Bacteria"],
            "phylum": ["Pseudomonadati"],
        }
    )

    result = oc.add_cellularity_smart(df)
    assert result.loc[0, "type"] == "Unicellular organism"


def test_add_cellularity_smart_distinguishes_eukaryotes() -> None:
    df = pd.DataFrame(
        {
            "genus": ["Homo", "Plasmodium"],
            "superkingdom": ["Eukaryota", "Eukaryota"],
            "phylum": ["Chordata", "Sar"],
        }
    )

    result = oc.add_cellularity_smart(df)
    assert result[result["genus"] == "Homo"].iloc[0]["type"] == "Multicellular organism"
    assert result[result["genus"] == "Plasmodium"].iloc[0]["type"] == "Unicellular organism"

