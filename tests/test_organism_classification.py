import pandas as pd

from library.organism_classification import (
    TYPE_MULTICELLULAR,
    TYPE_UNICELLULAR,
    TYPE_VIRAL,
    add_cellularity,
    add_cellularity_smart,
    classify_by_lineage,
    classify_record,
    normalize,
)


def test_normalize_handles_missing_values():
    assert normalize(None) == ""
    assert normalize("  ") == ""
    assert normalize("NA") == ""
    assert normalize(float("nan")) == ""
    assert normalize(pd.NA) == ""


def test_classify_by_lineage_prioritises_superkingdom():
    result = classify_by_lineage(
        genus="candida",
        superkingdom="viruses",
        phylum="",
        klass="",
    )
    assert result == TYPE_VIRAL


def test_classify_record_uses_all_lineage_levels():
    record = {
        "genus": "Chlamydomonas",
        "lineage_superkingdom": "Eukaryota",
        "lineage_phylum": "",
        "lineage_class": "",
    }
    assert classify_record(record) == TYPE_UNICELLULAR


def test_classify_by_phylum_and_class():
    assert (
        classify_by_lineage(
            genus=normalize(""),
            superkingdom=normalize("Eukaryota"),
            phylum=normalize("Apicomplexa"),
            klass=normalize(""),
        )
        == TYPE_UNICELLULAR
    )
    assert (
        classify_by_lineage(
            genus=normalize(""),
            superkingdom=normalize("Eukaryota"),
            phylum=normalize(""),
            klass=normalize("Saccharomycetes"),
        )
        == TYPE_UNICELLULAR
    )


def test_classify_returns_multicellular_by_default():
    assert (
        classify_by_lineage(
            genus="Homo",
            superkingdom="Eukaryota".lower(),
            phylum="Chordata",
            klass="Mammalia",
        )
        == TYPE_MULTICELLULAR
    )


def test_add_cellularity_adds_column_with_string_dtype():
    df = pd.DataFrame(
        {
            "genus": ["Candida", None],
            "lineage_superkingdom": ["Eukaryota", "Bacteria"],
            "lineage_phylum": ["", ""],
            "lineage_class": ["", ""],
        }
    )

    result = add_cellularity(df)
    assert result is not df
    assert list(result["type"]) == [TYPE_UNICELLULAR, TYPE_UNICELLULAR]
    assert str(result["type"].dtype) == "string"


def test_add_cellularity_smart_allows_custom_columns():
    df = pd.DataFrame(
        {
            "GenusName": ["Homo"],
            "SuperKingdom": ["Eukaryota"],
            "PhylumName": ["Chordata"],
            "ClassName": ["Mammalia"],
        }
    )

    result = add_cellularity_smart(
        df,
        genus_col="GenusName",
        superkingdom_col="SuperKingdom",
        phylum_col="PhylumName",
        class_col="ClassName",
        output_col="OrganismType",
    )

    assert list(result["OrganismType"]) == [TYPE_MULTICELLULAR]
    assert str(result["OrganismType"].dtype) == "string"
