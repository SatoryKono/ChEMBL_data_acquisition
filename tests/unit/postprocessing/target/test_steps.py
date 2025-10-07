from __future__ import annotations

import pandas as pd
import pytest

from library.postprocess.targets import steps
from library.postprocess.targets.steps import finalize_target_records, normalize_target_fields


@pytest.mark.unit
def test_normalize_target_fields__applies_taxonomy_and_identifier_normalization() -> None:
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
def test_normalize_target_fields__without_optional_flags_preserves_extra_columns() -> None:
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

    assert list(result.columns) == [
        "target_chembl_id",
        "pref_name",
        "target_type",
    ]
    assert pd.isna(result.loc[0, "target_type"])
    assert str(result["target_type"].dtype) == "string"
    assert str(result["pref_name"].dtype) == "string"


@pytest.mark.unit
def test_finalize_target_records__handles_empty_input_frame() -> None:
    frame = pd.DataFrame({"pref_name": []})

    result = finalize_target_records(frame)

    assert result.empty
    assert list(result.columns) == [
        "target_chembl_id",
        "pref_name",
        "target_type",
    ]
    assert str(result["target_chembl_id"].dtype) == "string"
    assert str(result["pref_name"].dtype) == "string"
    assert str(result["target_type"].dtype) == "string"


@pytest.mark.unit
def test_finalize_target_records__supports_sort_override() -> None:
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL2", "CHEMBL1"],
            "pref_name": ["Beta", "Alpha"],
            "target_type": ["type-b", "type-a"],
        }
    )

    result = finalize_target_records(frame, enforce_schema=False, sort_by=["pref_name"])

    assert result.loc[0, "pref_name"] == "Alpha"
    assert result.loc[1, "pref_name"] == "Beta"
    assert list(result.columns[:3]) == [
        "target_chembl_id",
        "pref_name",
        "target_type",
    ]


@pytest.mark.unit
def test_finalize_target_records__skips_schema_validation_when_disabled(monkeypatch) -> None:
    calls: list[pd.DataFrame] = []

    def _fake_validate(df: pd.DataFrame, *, context: str):  # pragma: no cover - sanity
        calls.append(df)
        return df

    monkeypatch.setattr(steps, "validate_targets", _fake_validate)

    frame = pd.DataFrame({"target_chembl_id": ["CHEMBL1"], "pref_name": ["Alpha"]})

    result = finalize_target_records(frame, enforce_schema=False)

    assert not calls
    assert list(result.columns)[:3] == [
        "target_chembl_id",
        "pref_name",
        "target_type",
    ]
