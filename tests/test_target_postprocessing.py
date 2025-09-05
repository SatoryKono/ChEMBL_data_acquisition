"""Tests for :mod:`library.target_postprocessing`.

The suite covers :func:`library.target_postprocessing.postprocess_targets` and
its convenience wrapper :func:`library.target_postprocessing.postprocess_file`.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library import target_postprocessing as tp


def test_postprocess_targets_merges_and_normalises() -> None:
    """``postprocess_targets`` merges fields and normalises tokens."""

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


def test_postprocess_file_roundtrip(tmp_path: Path) -> None:
    """Ensure ``postprocess_file`` writes the same result as ``postprocess_targets``."""

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

    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False)
    output_path = tmp_path / "out.csv"

    tp.postprocess_file(input_path, output_path)

    expected = tp.postprocess_targets(df).astype(str)
    result = pd.read_csv(output_path, dtype=str)
    pd.testing.assert_frame_equal(result, expected)
