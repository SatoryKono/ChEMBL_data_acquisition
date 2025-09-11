import pandas as pd
import pytest

from library.input_initialisation_library import (
    add_percentage,
    compute_status_statistics,
    get_percentage,
)


def test_get_percentage_basic() -> None:
    df = pd.DataFrame({"Filtered": ["good", "bad", "good"]})
    res = get_percentage(df, "activity")
    assert res.loc[res["Filtered"] == "good", "Count"].iat[0] == 2
    assert res.loc[res["Filtered"] == "Total", "Percentage, %"].iat[0] == 3


def test_get_percentage_missing_column() -> None:
    with pytest.raises(KeyError):
        get_percentage(pd.DataFrame({"foo": [1]}), "activity")


def test_add_percentage_and_compute_statistics() -> None:
    df = pd.DataFrame(
        {
            "activity_id": [1, 2, 3],
            "Filtered": ["good", "bad", "good"],
            "independent_IC50": [1, 0, 1],
            "non_independent_IC50": [0, 1, 0],
            "independent_Ki": [0, 0, 0],
            "non_independent_Ki": [0, 0, 0],
        }
    )
    res = compute_status_statistics(df, "activity")
    assert "activity.Percentage, %" in res.columns
    assert res.loc[res["Filtered"] == "good", "Count"].iat[0] == 2
    assert res.loc[res["Filtered"] == "Total", "Count"].iat[0] == 3


def test_add_percentage_renames_column() -> None:
    statistics = pd.DataFrame({"Filtered": ["good", "Total"], "metric": [1, 1]})
    percent_df = pd.DataFrame(
        {
            "Filtered": ["good", "Total"],
            "Count": [1, 2],
            "Percentage, %": [50.0, 100.0],
        }
    )
    res = add_percentage(statistics, percent_df, "activity")
    assert "activity.Percentage, %" in res.columns
    assert res.loc[res["Filtered"] == "good", "metric"].iat[0] == 1
