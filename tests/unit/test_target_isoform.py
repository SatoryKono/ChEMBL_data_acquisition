"""Unit tests for the isoform post-processing helpers."""

from __future__ import annotations

import pandas as pd
import pytest

from library.postprocessing.target import isoform


@pytest.mark.unit
def test_transform__supports_alias_columns():
    """Isoform transform accepts common legacy aliases for required columns."""

    frame = pd.DataFrame(
        {
            "isoform_synonyms": ["SynA|SynB"],
            "isoform_names": ["Isoform 1|Isoform 2"],
            "isoform_ids": ["P12345-1|P12345-2"],
            "primary_accession": ["P12345"],
            "chembl_id": ["CHEMBL123"],
        }
    )

    result = isoform._transform(frame).result

    expected = pd.DataFrame(
        [
            {
                "id": "P12345-1",
                "uniprot_id_primary": "P12345",
                "target_chembl_id": "CHEMBL123",
                "name": "isoform 1",
            },
            {
                "id": "P12345-1",
                "uniprot_id_primary": "P12345",
                "target_chembl_id": "CHEMBL123",
                "name": "syna",
            },
            {
                "id": "P12345-2",
                "uniprot_id_primary": "P12345",
                "target_chembl_id": "CHEMBL123",
                "name": "isoform 2",
            },
            {
                "id": "P12345-2",
                "uniprot_id_primary": "P12345",
                "target_chembl_id": "CHEMBL123",
                "name": "synb",
            },
        ]
    )

    pd.testing.assert_frame_equal(result, expected)


@pytest.mark.unit
def test_transform__missing_identifier_columns_warns_and_returns_empty_result():
    """A warning is emitted and an empty result returned when IDs are absent."""

    frame = pd.DataFrame(
        {
            "isoform_synonyms": ["SynA"],
            "isoform_names": ["Isoform 1"],
            "isoform_ids": ["P12345-1"],
        }
    )

    with pytest.warns(UserWarning) as captured:
        result = isoform._transform(frame)

    assert len(captured) == 1
    message = str(captured[0].message)
    assert "uniprot_id_primary" in message
    assert "target_chembl_id" in message

    for attribute in ("result", "combined", "dedup_stage1", "sorted_stage", "dedup_stage2"):
        frame_result = getattr(result, attribute)
        assert isinstance(frame_result, pd.DataFrame)
        assert frame_result.empty


@pytest.mark.unit
@pytest.mark.parametrize(
    "filename",
    [
        "out.csv",
        "out_chembl.csv",
        "out_normalized.csv",
        "out_uniprot.csv",
    ],
)
def test_matches_expected_input_name__supports_modern_exports(filename: str) -> None:
    """Modern ``out*.csv`` exports are accepted for isoform processing."""

    assert isoform._matches_expected_input_name(filename)
