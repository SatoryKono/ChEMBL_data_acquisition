import pandas as pd
import pytest

from library.schemas.normalize import normalize_testitems


@pytest.mark.unit
def test_normalize_testitems__standardises_relations_and_units() -> None:
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": [" CHEMBL1 ", "CHEMBL2"],
            "relation": ["<", ">"],
            "units": ["5 μM", "10 μM"],
        }
    )

    normalised = normalize_testitems(frame)

    assert list(normalised["molecule_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
    assert normalised["molecule_chembl_id"].dtype == "string"
    assert list(normalised["relation"]) == ["<=", ">="]
    assert list(normalised["units"]) == ["5 uM", "10 uM"]


@pytest.mark.unit
def test_normalize_testitems__preserves_non_string_values() -> None:
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "natural_product": [True],
            "measurement": [1.23],
        }
    )

    normalised = normalize_testitems(frame)

    assert normalised.equals(
        frame.assign(molecule_chembl_id=pd.Series(["CHEMBL1"], dtype="string"))
    )
