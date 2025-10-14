"""Smoke tests for activity output column filtering."""

from __future__ import annotations

import pandas as pd
import pytest
from scripts import get_activity_data


@pytest.mark.unit
@pytest.mark.pipeline_scenario("assembly")
def test_filter_activity_output_columns__drops_unwanted_columns_preserves_order():
    # Arrange
    source_columns = [
        "activity_id",
        "approx_cited_activity",
        "standard_value",
        "shuffled_cit",
        "target_chembl_id",
    ]
    frame = pd.DataFrame(
        {
            "activity_id": ["A1", "A2"],
            "approx_cited_activity": [True, False],
            "standard_value": [1.0, 2.0],
            "shuffled_cit": [pd.NA, pd.NA],
            "target_chembl_id": ["T1", "T2"],
        }
    )

    # Act
    filtered, dropped = get_activity_data._filter_activity_output_columns(
        frame,
        column_order=source_columns,
    )

    # Assert
    assert dropped == [
        "approx_cited_activity",
        "shuffled_cit",
    ]
    assert list(filtered.columns) == [
        "activity_id",
        "standard_value",
        "target_chembl_id",
    ]
    pd.testing.assert_frame_equal(
        filtered,
        frame[["activity_id", "standard_value", "target_chembl_id"]],
        check_dtype=True,
    )


@pytest.mark.unit
def test_filter_activity_output_columns__no_op_when_columns_absent():
    # Arrange
    frame = pd.DataFrame(
        {
            "activity_id": ["A1"],
            "standard_value": [1.0],
        }
    )

    # Act
    filtered, dropped = get_activity_data._filter_activity_output_columns(
        frame,
        column_order=["activity_id", "standard_value"],
    )

    # Assert
    assert dropped == []
    assert list(filtered.columns) == ["activity_id", "standard_value"]
    pd.testing.assert_frame_equal(filtered, frame, check_dtype=True)
