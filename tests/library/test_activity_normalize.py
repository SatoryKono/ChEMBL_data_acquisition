"""Unit tests for :mod:`library.postprocessing.activity.steps`."""

from __future__ import annotations

import pandas as pd
from pandas.testing import assert_frame_equal

from library.postprocessing.activity.steps import (
    ACTIVITY_COLUMN_ORDER,
    normalize_activity_frame,
)


def test_normalize_activity_frame_enforces_order_and_types() -> None:
    raw = pd.DataFrame(
        [
            {
                "activity_id": " 2 ",
                "assay_chembl_id": " CHEMBL2 ",
                "target_chembl_id": "",
                "standard_type": " ic50 ",
                "standard_relation": " = ",
                "standard_value": "42.5",
                "standard_units": "nm",
            },
            {
                "activity_id": "1",
                "assay_chembl_id": None,
                "target_chembl_id": "CHEMBL3",
                "standard_type": "Ki",
                "standard_relation": "",
                "standard_value": "not-a-number",
                "standard_units": None,
            },
        ]
    )

    normalised = normalize_activity_frame(raw)

    assert tuple(normalised.columns) == ACTIVITY_COLUMN_ORDER
    assert list(normalised["activity_id"]) == [1, 2]
    assert normalised["activity_id"].dtype == "Int64"
    assert normalised["standard_value"].dtype == "Float64"
    assert normalised["standard_units"].dtype == pd.StringDtype()
    assert pd.isna(normalised.loc[0, "standard_units"])
    assert normalised.loc[1, "standard_units"] == "NM"
    assert normalised.loc[1, "standard_relation"] == "="
    assert pd.isna(normalised.loc[0, "standard_relation"])
    assert pd.isna(normalised.loc[0, "standard_value"])


def test_normalize_activity_frame_handles_missing_columns() -> None:
    raw = pd.DataFrame({"activity_id": ["5"], "standard_units": [" uM "]})

    normalised = normalize_activity_frame(raw)

    expected = pd.DataFrame(
        {
            "activity_id": pd.Series([5], dtype="Int64"),
            "assay_chembl_id": pd.Series([pd.NA], dtype=pd.StringDtype()),
            "target_chembl_id": pd.Series([pd.NA], dtype=pd.StringDtype()),
            "standard_type": pd.Series([pd.NA], dtype=pd.StringDtype()),
            "standard_relation": pd.Series([pd.NA], dtype=pd.StringDtype()),
            "standard_value": pd.Series([pd.NA], dtype="Float64"),
            "standard_units": pd.Series(["UM"], dtype=pd.StringDtype()),
        }
    )

    assert_frame_equal(normalised, expected)

