from __future__ import annotations

import pandas as pd
import pytest

from library.schemas.normalize import normalize_testitems


@pytest.mark.unit
@pytest.mark.parametrize(
    "relation, expected",
    [("<", "<="), (">", ">="), ("<=", "<="), ("invalid", "invalid")],
)
def test_normalize_testitems__relation_mapping(relation: str, expected: str) -> None:
    frame = pd.DataFrame({"relation": [relation]})
    normalised = normalize_testitems(frame)
    assert normalised.loc[0, "relation"] == expected


@pytest.mark.unit
@pytest.mark.pipeline_scenario("normalization")
def test_normalize_testitems__identifier_columns_trimmed() -> None:
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": [" cheMBL1 "],
            "parent_molecule_chembl_id": ["  chembl2"],
        }
    )
    normalised = normalize_testitems(frame)
    assert normalised.loc[0, "molecule_chembl_id"] == "cheMBL1"
    assert normalised["molecule_chembl_id"].dtype == "string"
    assert normalised.loc[0, "parent_molecule_chembl_id"] == "chembl2"


@pytest.mark.unit
def test_normalize_testitems__micro_sign_replaced() -> None:
    frame = pd.DataFrame({"units": ["5 μM", "10 uM"]})
    normalised = normalize_testitems(frame)
    assert list(normalised["units"]) == ["5 uM", "10 uM"]


@pytest.mark.unit
def test_normalize_testitems__preserves_boolean_dtype() -> None:
    frame = pd.DataFrame({"natural_product": pd.Series([True, pd.NA], dtype="boolean")})
    normalised = normalize_testitems(frame)
    assert normalised["natural_product"].dtype == "boolean"
    assert normalised["natural_product"].tolist() == [True, pd.NA]


@pytest.mark.unit
@pytest.mark.pipeline_scenario("missing_data")
def test_normalize_testitems__retains_missing_values() -> None:
    frame = pd.DataFrame({"molecule_chembl_id": [pd.NA]})
    normalised = normalize_testitems(frame)
    assert pd.isna(normalised.loc[0, "molecule_chembl_id"])


@pytest.mark.unit
def test_normalize_testitems__does_not_mutate_original() -> None:
    frame = pd.DataFrame({"relation": ["<"]})
    _ = normalize_testitems(frame)
    assert frame.loc[0, "relation"] == "<"
