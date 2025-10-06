"""Smoke tests for ensuring unwanted assay columns are excluded from output."""

from __future__ import annotations

import pandas as pd
import pytest

from scripts.get_assay_data import ASSAY_OUTPUT_DROP_COLUMNS, remove_assay_output_columns


@pytest.mark.unit
def test_remove_assay_output_columns__drops_specified_fields() -> None:
    """The helper removes only the disallowed assay columns preserving order."""

    columns = [
        "assay_chembl_id",
        "ASSAY_ID",
        "Target TYPE",
        "assay_type",
        "acts_per_assay_step5",
        "custom_field",
    ]
    values = list(range(len(columns)))
    frame = pd.DataFrame([values], columns=columns)

    result = remove_assay_output_columns(frame)

    expected_columns = [
        column for column in columns if column not in ASSAY_OUTPUT_DROP_COLUMNS
    ]

    assert list(result.columns) == expected_columns
    pd.testing.assert_frame_equal(
        result,
        pd.DataFrame([[0, 3, 5]], columns=expected_columns),
        check_dtype=False,
    )
    assert list(frame.columns) == columns, "original frame must remain unchanged"


@pytest.mark.unit
def test_remove_assay_output_columns__handles_missing_fields() -> None:
    """Columns absent from the dataframe are ignored without raising errors."""

    frame = pd.DataFrame(
        [
            {
                "assay_chembl_id": "CHEMBL1",
                "assay_type": "binding",
                "custom_field": "value",
            }
        ]
    )

    result = remove_assay_output_columns(frame)

    assert list(result.columns) == ["assay_chembl_id", "assay_type", "custom_field"]
    pd.testing.assert_frame_equal(result, frame, check_like=True)
