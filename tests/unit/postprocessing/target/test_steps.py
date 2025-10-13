from __future__ import annotations

import json

import pandas as pd
import pytest

from library.postprocessing.targets.steps import (
    enrich_target_synonyms,
    finalize_target_records,
    normalize_target_fields,
)


@pytest.mark.unit
def test_normalize_target_fields__applies_taxonomy_and_identifier_normalization() -> (
    None
):
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["  chembl123  ", None],
            "pref_name": ["  Alpha  ", "Beta"],
            "organism": ["  Homo sapiens  ", "  "],
            "lineage_superkingdom": ["  Eukaryota  ", None],
            "lineage_phylum": [None, "  Chordata"],
            "taxon_id": [" 9606 ", ""],
            "tax_id": ["9606", "  "],
            "target_id": ["  chembl123  ", "CHEMBL999"],
        }
    )

    result = normalize_target_fields(
        frame,
        normalize_taxonomy=True,
        fill_missing_identifiers=True,
    )

    assert result.loc[0, "target_chembl_id"] == "CHEMBL123"
    assert result.loc[1, "target_chembl_id"] == "CHEMBL999"
    assert result.loc[0, "organism"] == "Homo sapiens"
    assert pd.isna(result.loc[1, "organism"])
    assert result.loc[0, "lineage_superkingdom"] == "Eukaryota"
    assert result.loc[1, "lineage_phylum"] == "Chordata"
    assert result["taxon_id"].tolist() == [9606, pd.NA]
    assert result["tax_id"].tolist() == [9606, pd.NA]

    # ensure the original dataframe was not mutated
    assert frame.loc[0, "target_chembl_id"] == "  chembl123  "
    assert frame.loc[1, "target_chembl_id"] is None


@pytest.mark.unit
def test_normalize_target_fields__derives_classifications_and_synonyms() -> None:
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "protein_classifications": [
                json.dumps(
                    [
                        {
                            "protein_classification": {
                                "pref_name": "Enzyme",
                                "class_level": 1,
                            }
                        },
                        {
                            "protein_classification": {
                                "pref_name": "Kinase family",
                                "class_level": 2,
                            }
                        },
                    ]
                ),
                "",
            ],
            "protein_class_pred_L1": [None, "Other Protein Target"],
            "protein_class_pred_L2": [None, "Miscellaneous"],
            "protein_synonym_list": ["Alpha|Beta", None],
            "protein_name_canonical": ["Canonical A", "Canonical B"],
            "protein_name_alt": ["AltA", None],
            "pref_name": [" Alpha ", "Beta"],
            "gene_symbol": ["GENA", None],
            "gtop_synonyms": ["delta;epsilon", None],
        }
    )

    result = normalize_target_fields(frame)

    assert result["target_class"].tolist() == ["Enzyme", "Other Protein Target"]
    assert result["protein_family"].tolist() == ["Kinase family", "Miscellaneous"]

    chembl_tokens = result.loc[0, "chembl_synonyms"].split("|")
    assert chembl_tokens == ["Alpha", "Beta", "Canonical A", "AltA", "GENA"]
    assert result.loc[0, "gtopdb_synonyms"] == "delta|epsilon"

    assert result.loc[1, "chembl_synonyms"] == "Canonical B|Beta"
    assert pd.isna(result.loc[1, "gtopdb_synonyms"])


@pytest.mark.unit
def test_normalize_target_fields__without_optional_flags_preserves_extra_columns() -> (
    None
):
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["chembl1"],
            "pref_name": ["  Gamma"],
            "taxon_id": [" 9606 "],
        }
    )

    result = normalize_target_fields(frame)

    assert result.loc[0, "target_chembl_id"] == "CHEMBL1"
    assert result.loc[0, "pref_name"] == "Gamma"
    assert result.loc[0, "taxon_id"] == " 9606 "


@pytest.mark.unit
def test_finalize_target_records__fills_missing_required_columns() -> None:
    frame = pd.DataFrame({"target_chembl_id": ["CHEMBL1"], "pref_name": ["Alpha"]})

    result = finalize_target_records(frame)

    assert list(result.columns[:3]) == [
        "target_chembl_id",
        "pref_name",
        "target_type",
    ]
    for column in [
        "organism",
        "target_class",
        "protein_family",
        "synonyms",
        "pipeline_version",
    ]:
        assert column in result.columns
    assert pd.isna(result.loc[0, "target_type"])
    assert str(result["target_type"].dtype) == "string"
    assert str(result["pref_name"].dtype) == "string"


@pytest.mark.unit
def test_finalize_target_records__handles_empty_input_frame() -> None:
    frame = pd.DataFrame({"pref_name": []})

    result = finalize_target_records(frame)

    assert result.empty
    assert list(result.columns[:3]) == [
        "target_chembl_id",
        "pref_name",
        "target_type",
    ]
    for column in [
        "organism",
        "target_class",
        "protein_family",
        "synonyms",
        "pipeline_version",
    ]:
        assert column in result.columns
    assert str(result["target_chembl_id"].dtype) == "string"
    assert str(result["pref_name"].dtype) == "string"
    assert str(result["target_type"].dtype) == "string"


@pytest.mark.unit
def test_enrich_target_synonyms__combines_sources_with_separator() -> None:
    frame = pd.DataFrame(
        {
            "synonyms": ["beta | Alpha", None],
            "chembl_synonyms": ["alpha; delta", "theta"],
            "gtopdb_synonyms": ["gamma, beta", ""],
        }
    )

    result = enrich_target_synonyms(
        frame,
        synonym_sources=["chembl", "gtopdb"],
        preferred_separator="; ",
    )

    assert result["synonyms"].tolist() == ["beta; Alpha; delta; gamma", "theta"]


@pytest.mark.unit
def test_enrich_target_synonyms__ignores_unrecognised_parameters() -> None:
    frame = pd.DataFrame(
        {
            "synonyms": ["alpha"],
            "chembl_synonyms": ["beta"],
        }
    )

    result = enrich_target_synonyms(
        frame,
        synonym_sources=["chembl"],
        preferred_separator="; ",
        future_flag=True,
    )

    assert result["synonyms"].tolist() == ["alpha; beta"]


@pytest.mark.unit
def test_finalize_target_records__supports_optional_flags(monkeypatch) -> None:
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL2", "CHEMBL1"],
            "pref_name": ["Beta", "Alpha"],
            "target_type": ["protein", "protein"],
        }
    )

    monkeypatch.setattr(
        "library.postprocessing.targets.steps.validate_targets",
        lambda *args, **kwargs: pytest.fail("validate_targets should not be called"),
    )

    result = finalize_target_records(
        frame,
        enforce_schema=False,
        sort_by=["target_chembl_id"],
    )

    assert result["target_chembl_id"].tolist() == ["CHEMBL1", "CHEMBL2"]


@pytest.mark.unit
def test_finalize_target_records__populates_target_type_from_relationship() -> None:
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL42"],
            "pref_name": ["Some target"],
            "relationship": ["SINGLE PROTEIN"],
        }
    )

    result = finalize_target_records(frame)

    assert result.loc[0, "target_type"] == "SINGLE PROTEIN"


@pytest.mark.unit
def test_finalize_target_records__fills_target_type_from_dictionary(
    monkeypatch,
) -> None:
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["chembl1", "CHEMBL2"],
            "pref_name": ["Alpha", "Beta"],
            "target_type": [pd.NA, None],
        }
    )

    lookup = pd.Series(
        {
            "CHEMBL1": "Multicellular organism",
            "CHEMBL2": "Unicellular organism",
        },
        dtype="string",
    )

    monkeypatch.setattr(
        "library.postprocessing.targets.steps._get_target_type_lookup",
        lambda: lookup,
    )

    result = finalize_target_records(frame)

    mapping = {
        key.upper(): value
        for key, value in result.set_index("target_chembl_id")["target_type"].items()
    }

    assert mapping["CHEMBL1"] == "Multicellular organism"
    assert mapping["CHEMBL2"] == "Unicellular organism"


def test_finalize_target_records__populates_target_type_from_description() -> None:
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL9000"],
            "pref_name": ["Example target"],
            "target_type_description": ["PROTEIN COMPLEX"],
        }
    )

    result = finalize_target_records(frame)

    assert result.loc[0, "target_type"] == "PROTEIN COMPLEX"
    assert str(result["target_type"].dtype) == "string"
