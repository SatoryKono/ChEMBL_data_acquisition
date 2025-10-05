"""Unit tests for the isoform post-processing helpers."""

from __future__ import annotations

from pathlib import Path

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
def test_transform__missing_identifier_columns_raises_key_error():
    """A descriptive error is raised when identifier columns are absent."""

    frame = pd.DataFrame(
        {
            "isoform_synonyms": ["SynA"],
            "isoform_names": ["Isoform 1"],
            "isoform_ids": ["P12345-1"],
        }
    )

    with pytest.raises(KeyError) as exc:
        isoform._transform(frame)

    message = str(exc.value)
    assert "uniprot_id_primary" in message
    assert "target_chembl_id" in message


@pytest.mark.unit
def test_matches_expected_input_name__accepts_custom_normalized() -> None:
    """Custom stems ending with ``_normalized`` are accepted."""

    assert isoform._matches_expected_input_name("custom_output_normalized.csv")
    assert not isoform._matches_expected_input_name("custom_output_chembl.csv")


@pytest.mark.unit
def test_process_targets__custom_normalized_name(tmp_path: Path, snapshot_resource: Path) -> None:
    """``process_targets`` handles user-specified ``*_normalized.csv`` exports."""

    source = snapshot_resource / "target_isoform_minimal.csv"
    input_path = tmp_path / "user_defined_normalized.csv"
    input_path.write_bytes(source.read_bytes())

    output_path = isoform.process_targets(str(input_path), verbose=False)

    assert output_path.name == "isoform.user_defined_normalized.csv"
    assert output_path.parent == input_path.parent
