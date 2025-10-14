from __future__ import annotations

import pandas as pd

from library.postprocessing.target.steps import enrich_target_metadata


def test_enrich_target_metadata_respects_existing_values() -> None:
    base = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "pref_name": ["Target 1", "Target 2"],
            "target_type": ["protein", "protein"],
            "organism": ["Human", "Human"],
            "uniprot_id": ["P12345", "P67890"],
            "gene_symbol": ["GENE1", "GENE2"],
            "target_class": [pd.NA, "Existing"],
            "protein_family": [pd.NA, "FamilyB"],
            "synonyms": [pd.NA, "BaseSyn"],
        }
    ).astype("string")

    uniprot_metadata = pd.DataFrame(
        {
            "uniprot_id": ["P12345", "P67890"],
            "protein_family": ["Kinase", "Override"],
            "synonyms": ["MetaSyn1", "MetaSyn2"],
        }
    ).astype("string")

    gtopdb_metadata = pd.DataFrame(
        {
            "uniprot_id": ["P12345", "P67890"],
            "target_class": ["ClassA", "ClassB"],
        }
    ).astype("string")

    result = enrich_target_metadata(base, uniprot_metadata, gtopdb_metadata)

    assert result.loc[0, "target_class"] == "ClassA"
    assert result.loc[1, "target_class"] == "Existing"
    assert result.loc[0, "protein_family"] == "Kinase"
    assert result.loc[1, "protein_family"] == "FamilyB"
    assert result.loc[0, "synonyms"] == "MetaSyn1"
    assert result.loc[1, "synonyms"] == "BaseSyn"
    assert result["uniprot_id"].dtype == pd.StringDtype()
