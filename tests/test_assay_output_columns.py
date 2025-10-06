"""Smoke tests for assay output column filtering."""

from __future__ import annotations

import logging

import pandas as pd
import pytest

from scripts.get_assay_data import _drop_assay_output_columns


@pytest.mark.unit
def test_drop_assay_output_columns__removes_and_preserves_order(caplog: pytest.LogCaptureFixture) -> None:
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
    frame = pd.DataFrame([
        ["A", "legacy-1", "type-1", "primary", "1.0", "substrate", "value"],
        ["B", "legacy-2", "type-2", "backup", "2.0", "substrate", "other"],
    ], columns=source_columns)

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
