"""Tests for :mod:`library.target_postprocessing`."""

from __future__ import annotations

import pandas as pd

from library import target_postprocessing as tp


def test_postprocess_targets_merges_and_normalises() -> None:
    df = pd.DataFrame(
        {
            "chembl_id": ["CHEMBL1"],
            "uniProtkbId": ["P12345-2_SOMETHING"],
            "uniprot_id": ["P12345"],
            "secondaryAccessions": ["S12345|S67890"],
            "pref_name": ["Protein kinase"],
            "gene_name_x": ["g1"],
            "geneName": [""],
            "gene": ["g2|g3"],
            "secondaryAccessionNames": ["name2|name3"],
            "component_description": ["desc"],
            "chembl_alternative_name": ["alt"],
            "recommendedName": ["Rec Name"],
            "names": ["name1|name4"],
            "synonyms_x": ["synX"],
            "synonyms": ["syn1|syn2"],
            "ec_number": ["1.1.1.1"],
            "ec_code": ["2.2.2.2"],
        }
    )
    out = tp.postprocess_targets(df)

    row = out.iloc[0]
    assert row["uniprotkb_Id"] == "P12345"
    assert row["secondary_uniprot_id"] == "S12345|S67890"
    assert row["recommended_name"] == "Protein kinase"
    assert row["gene_name"] == "G1"
    assert row["ec_number"] == "1.1.1.1|2.2.2.2"
    assert row["isoform_names"] == "-"
    assert row["synonyms"].startswith("g2|g3|g1|name2|name3")
