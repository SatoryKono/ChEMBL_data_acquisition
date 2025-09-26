"""Tests for :mod:`library.organism_classification`."""

from __future__ import annotations

import pandas as pd

from library.organism_classification import (
    CELLULARITY_MULTICELLULAR,
    CELLULARITY_UNICELLULAR,
    CELLULARITY_VIRUS,
    OrganismClassificationRules,
    add_cellularity,
    add_cellularity_smart,
    classify_by_lineage,
    classify_record,
    normalize,
)


def test_normalize_handles_missing_tokens() -> None:
    assert normalize("  Homo  ") == "Homo"
    assert normalize("") is None
    assert normalize(None) is None
    assert normalize("NaN") is None
    assert normalize(pd.NA) is None


def test_classify_by_lineage_priority() -> None:
    rules = OrganismClassificationRules()

    assert (
        classify_by_lineage("Viruses", "Riboviria", None, rules=rules)
        == CELLULARITY_VIRUS
    )
    assert (
        classify_by_lineage("Eukaryota", "Viridiplantae", "Chlorophyceae", rules=rules)
        == CELLULARITY_UNICELLULAR
    )
    assert (
        classify_by_lineage("Eukaryota", "Chordata", "Mammalia", rules=rules)
        == CELLULARITY_MULTICELLULAR
    )
    assert classify_by_lineage(None, None, None, rules=rules) is None


def test_classify_record_respects_custom_columns() -> None:
    record = {"super": "Bacteria", "phy": "Proteobacteria", "cls": "Gammaproteobacteria"}
    assert (
        classify_record(
            record,
            superkingdom_col="super",
            phylum_col="phy",
            class_col="cls",
        )
        == CELLULARITY_UNICELLULAR
    )


def test_add_cellularity_appends_string_column() -> None:
    df = pd.DataFrame(
        {
            "lineage_superkingdom": ["Eukaryota", "Viruses", None],
            "lineage_phylum": ["Chordata", "Riboviria", "Sar"],
            "lineage_class": ["Mammalia", None, None],
        }
    )

    result = add_cellularity(df)

    assert result is not df
    assert result["cellularity"].dtype == "string"
    assert list(result["cellularity"].astype(object)) == [
        CELLULARITY_MULTICELLULAR,
        CELLULARITY_VIRUS,
        CELLULARITY_UNICELLULAR,
    ]


def test_add_cellularity_smart_preserves_existing_labels() -> None:
    df = pd.DataFrame(
        {
            "lineage_superkingdom": ["Eukaryota", "Bacteria", None],
            "lineage_phylum": ["Viridiplantae", "Proteobacteria", None],
            "lineage_class": ["Chlorophyceae", None, None],
            "organism_type": ["Custom", "", pd.NA],
        }
    )

    result = add_cellularity_smart(df)

    values = result["organism_type"].astype("string")
    assert list(values.iloc[:2].astype(object)) == [
        "Custom",
        CELLULARITY_UNICELLULAR,
    ]
    assert values.iloc[2] is pd.NA
