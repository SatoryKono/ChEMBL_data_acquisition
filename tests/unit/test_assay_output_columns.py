"""Smoke tests for ensuring unwanted assay columns are excluded from output."""

from __future__ import annotations

import logging

import pandas as pd
import pytest

from library.schemas.assays import AssaysSchema
from scripts.get_assay_data import (
    ASSAY_OUTPUT_DROP_COLUMNS,
    _drop_assay_output_columns,
    remove_assay_output_columns,
)


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


@pytest.mark.unit
def test_drop_assay_output_columns__removes_and_preserves_order(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """The whitelist keeps column order while dropping the disallowed set."""

    source_columns = [
        "assay_chembl_id",
        "ASSAY_ID",
        "assay_type",
        "Target TYPE",
        "version",
        "substrate_name",
        "custom_field",
    ]
    frame = pd.DataFrame(
        [
            ["A", "legacy-1", "type-1", "primary", "1.0", "substrate", "value"],
            ["B", "legacy-2", "type-2", "backup", "2.0", "substrate", "other"],
        ],
        columns=source_columns,
    )

    caplog.set_level(logging.INFO)

    filtered = _drop_assay_output_columns(frame)

    assert list(filtered.columns) == [
        "assay_chembl_id",
        "assay_type",
        "custom_field",
    ]
    pd.testing.assert_frame_equal(
        filtered,
        pd.DataFrame(
            [
                ["A", "type-1", "value"],
                ["B", "type-2", "other"],
            ],
            columns=["assay_chembl_id", "assay_type", "custom_field"],
        ),
    )

    messages = [record.getMessage() for record in caplog.records]
    assert any(
        "Dropped columns from output.assay_*: ASSAY_ID, Target TYPE, substrate_name, version"
        in message
        for message in messages
    )


@pytest.mark.unit
def test_drop_assay_output_columns__matches_schema_columns() -> None:
    """The helper preserves exactly the columns defined in ``AssaysSchema``."""

    schema_columns = list(AssaysSchema.columns.keys())

    # Safety net: if a dropped column sneaks into the schema we want the test to fail
    # with a clear message rather than producing duplicate columns.
    overlap = set(schema_columns) & set(ASSAY_OUTPUT_DROP_COLUMNS)
    assert not overlap, f"Schema unexpectedly includes dropped columns: {sorted(overlap)}"

    input_columns = schema_columns + ASSAY_OUTPUT_DROP_COLUMNS
    sample_row = {column: f"value-{idx}" for idx, column in enumerate(input_columns)}
    frame = pd.DataFrame([sample_row], columns=input_columns)

    filtered = _drop_assay_output_columns(frame)

    assert list(filtered.columns) == schema_columns
