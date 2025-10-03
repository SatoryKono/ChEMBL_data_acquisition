import pandas as pd
import pytest

from library.pipelines.testitem import cli


@pytest.mark.unit
def test_integrate_missing_identifiers__appends_placeholders_in_order() -> None:
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL2"],
            "value": [10],
        }
    )
    missing = ["CHEMBL3"]
    requested = ["CHEMBL2", "CHEMBL1", "CHEMBL3"]

    result = cli.integrate_missing_identifiers(
        df,
        missing_ids=missing,
        requested_ids=requested,
    )

    assert list(result["molecule_chembl_id"]) == ["CHEMBL2", "CHEMBL3"]
    assert pd.isna(result.loc[1, "value"])


@pytest.mark.unit
def test_integrate_missing_identifiers__restores_requested_order() -> None:
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL3", "CHEMBL1"],
        }
    )
    missing: list[str] = []
    requested = ["CHEMBL1", "CHEMBL2", "CHEMBL3"]

    result = cli.integrate_missing_identifiers(
        df,
        missing_ids=missing,
        requested_ids=requested,
    )

    assert list(result["molecule_chembl_id"]) == ["CHEMBL1", "CHEMBL3"]


@pytest.mark.unit
def test_integrate_missing_identifiers__ignores_blank_requested_ids() -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    result = cli.integrate_missing_identifiers(
        df,
        missing_ids=["CHEMBL2"],
        requested_ids=["", None, "CHEMBL1", "CHEMBL2"],
    )

    assert list(result["molecule_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]


@pytest.mark.unit
def test_integrate_missing_identifiers__missing_ids_appended_at_end() -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL3"]})
    result = cli.integrate_missing_identifiers(
        df,
        missing_ids=["CHEMBL1"],
        requested_ids=["CHEMBL3"],
    )

    assert list(result["molecule_chembl_id"]) == ["CHEMBL3", "CHEMBL1"]
