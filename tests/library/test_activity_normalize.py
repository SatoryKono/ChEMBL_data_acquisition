"""Unit tests covering deterministic activity normalization helpers."""

from __future__ import annotations

import pandas as pd

from library.postprocessing.activity.steps import (
    ACTIVITY_COLUMN_ORDER,
    ACTIVITY_SORT_COLUMNS,
    build_activity_reports,
    normalize_activity_frame,
)


def test_normalize_activity_frame_enforces_schema_and_order() -> None:
    raw = pd.DataFrame(
        [
            {
                "activity_id": "2",
                "assay_chembl_id": " chembl2 ",
                "target_chembl_id": " ",
                "standard_type": " ic50 ",
                "standard_relation": " ~= ",
                "standard_value": "1.23",
                "standard_units": "nm",
            },
            {
                "activity_id": "1",
                "assay_chembl_id": "chembl1",
                "target_chembl_id": "CHEMBL19999",
                "standard_type": None,
                "standard_relation": "=",
                "standard_value": "",
                "standard_units": "um",
            },
        ]
    )

    normalized = normalize_activity_frame(raw)

    assert list(normalized.columns) == list(ACTIVITY_COLUMN_ORDER)
    assert tuple(normalized.index) == (0, 1)
    assert tuple(normalized["activity_id"]) == (1, 2)
    assert str(normalized["activity_id"].dtype) == "Int64"
    assert str(normalized["standard_value"].dtype) == "Float64"

    expected_assays = ("CHEMBL1", "CHEMBL2")
    assert tuple(normalized["assay_chembl_id"]) == expected_assays

    assert pd.isna(normalized.loc[1, "target_chembl_id"])
    assert normalized.loc[0, "target_chembl_id"] == "CHEMBL19999"
    assert normalized.loc[0, "standard_units"] == "UM"
    assert normalized.loc[1, "standard_units"] == "NM"
    assert normalized.loc[0, "standard_relation"] == "="
    assert normalized.loc[1, "standard_relation"] == "~="
    assert normalized.loc[1, "standard_type"] == "IC50"

    # Ensure deterministic ordering respects the configured sort columns.
    assert ACTIVITY_SORT_COLUMNS == (
        "activity_id",
        "assay_chembl_id",
        "target_chembl_id",
        "standard_type",
        "standard_relation",
    )


def test_build_activity_reports_handles_small_frames() -> None:
    frame = pd.DataFrame(
        [
            {
                "activity_id": 3,
                "assay_chembl_id": "CHEMBL100",
                "target_chembl_id": "CHEMBL200",
                "standard_type": "IC50",
                "standard_relation": "=",
                "standard_value": 5.0,
                "standard_units": "NM",
            },
            {
                "activity_id": 4,
                "assay_chembl_id": "CHEMBL101",
                "target_chembl_id": "CHEMBL201",
                "standard_type": "IC50",
                "standard_relation": "=",
                "standard_value": 7.0,
                "standard_units": "NM",
            },
        ],
        columns=ACTIVITY_COLUMN_ORDER,
    )

    quality, correlation = build_activity_reports(frame)

    assert "column" in quality.columns
    assert tuple(quality["column"]) == tuple(sorted(quality["column"]))
    if not correlation.empty:
        assert list(correlation.index) == list(correlation.columns)

    empty_quality, empty_correlation = build_activity_reports(pd.DataFrame(columns=ACTIVITY_COLUMN_ORDER))
    assert empty_quality.empty
    assert empty_correlation.empty
