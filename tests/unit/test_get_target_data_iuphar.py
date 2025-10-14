import pandas as pd
import pytest
from scripts.get_target_data import _prepare_iuphar_merge_frames


@pytest.mark.unit
def test_prepare_iuphar_merge_frames__overrides_iuphar_columns():
    combined = pd.DataFrame(
        {
            "uniprot_id": ["P1", "P2"],
            "target_id": ["", ""],
            "GuidetoPHARMACOLOGY": ["", ""],
            "other": ["x", "y"],
        }
    )
    iuphar = pd.DataFrame(
        {
            "uniprot_id": ["P1", "P2"],
            "target_id": ["T1", "T2"],
            "GuidetoPHARMACOLOGY": ["G1", "G2"],
            "IUPHAR_family_id": ["F1", "F2"],
            "IUPHAR_type": ["Type1", "Type2"],
            "extra": ["foo", "bar"],
        }
    )

    updated_combined, trimmed_iuphar = _prepare_iuphar_merge_frames(combined, iuphar)

    assert "target_id" not in updated_combined.columns
    assert "GuidetoPHARMACOLOGY" not in updated_combined.columns
    assert trimmed_iuphar.columns.tolist() == [
        "uniprot_id",
        "target_id",
        "GuidetoPHARMACOLOGY",
        "IUPHAR_family_id",
        "IUPHAR_type",
        "extra",
    ]
    assert trimmed_iuphar.loc[0, "target_id"] == "T1"
    assert trimmed_iuphar.loc[0, "IUPHAR_family_id"] == "F1"


@pytest.mark.unit
def test_prepare_iuphar_merge_frames__keeps_non_override_columns_only_once():
    combined = pd.DataFrame(
        {
            "uniprot_id": ["P1"],
            "shared": ["chembl"],
        }
    )
    iuphar = pd.DataFrame(
        {
            "uniprot_id": ["P1"],
            "shared": ["iuphar"],
            "IUPHAR_class": ["Class"],
        }
    )

    updated_combined, trimmed_iuphar = _prepare_iuphar_merge_frames(combined, iuphar)

    assert "shared" in updated_combined.columns
    assert trimmed_iuphar.columns.tolist() == ["uniprot_id", "IUPHAR_class"]
