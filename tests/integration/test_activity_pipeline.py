"""Integration tests for the streamlined activity pipeline."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from library.postprocessing.activity.steps import (
    ACTIVITY_COLUMN_ORDER,
    ActivityData,
    fetch_normalize_activity,
    generate_activity_reports,
)


class _StubChemblClient:
    """Stub Chembl client returning pre-defined activity pages."""

    def __init__(self, payloads: Sequence[dict[str, object]]) -> None:
        self.payloads = list(payloads)
        self.calls: list[dict[str, object]] = []

    def list_activities(
        self,
        *,
        limit: int,
        offset: int = 0,
        fields: Sequence[str] | None = None,
        params=None,
        timeout=None,
    ) -> dict[str, object]:
        call = {
            "limit": limit,
            "offset": offset,
            "fields": tuple(fields or ()),
        }
        self.calls.append(call)
        if self.payloads:
            return self.payloads.pop(0)
        return {"activities": []}


def test_fetch_normalize_activity_combines_pages_and_validates() -> None:
    payloads = [
        {
            "activities": [
                {
                    "activity_id": "3",
                    "assay_chembl_id": "A3",
                    "target_chembl_id": "T2",
                    "standard_type": "IC50",
                    "standard_relation": "=",
                    "standard_value": "10",
                    "standard_units": "nm",
                },
                {
                    "activity_id": "1",
                    "assay_chembl_id": "A1",
                    "target_chembl_id": "T1",
                    "standard_type": "EC50",
                    "standard_relation": ">",
                    "standard_value": "5",
                    "standard_units": "nm",
                },
            ],
            "page_meta": {"next": "next-page"},
        },
        {
            "activities": [
                {
                    "activity_id": "2",
                    "assay_chembl_id": "A2",
                    "target_chembl_id": "T3",
                    "standard_type": "Ki",
                    "standard_relation": "=",
                    "standard_value": "1.5",
                    "standard_units": "um",
                }
            ],
            "page_meta": {},
        },
    ]
    client = _StubChemblClient(payloads)

    frame = fetch_normalize_activity(3, chembl_client=client)

    assert list(frame["activity_id"]) == [1, 2, 3]
    assert tuple(frame.columns) == ACTIVITY_COLUMN_ORDER
    assert client.calls[0] == {"limit": 3, "offset": 0, "fields": _expected_fields()}
    assert client.calls[1] == {"limit": 1, "offset": 2, "fields": _expected_fields()}


def _expected_fields() -> tuple[str, ...]:
    return (
        "activity_id",
        "assay_chembl_id",
        "target_chembl_id",
        "standard_type",
        "standard_relation",
        "standard_value",
        "standard_units",
    )


def test_generate_activity_reports_produces_qc_artifacts() -> None:
    frame = pd.DataFrame(
        {
            "activity_id": pd.Series([1, 2], dtype="Int64"),
            "assay_chembl_id": pd.Series(["A1", "A2"], dtype=pd.StringDtype()),
            "target_chembl_id": pd.Series(["T1", "T2"], dtype=pd.StringDtype()),
            "standard_type": pd.Series(["IC50", "IC50"], dtype=pd.StringDtype()),
            "standard_relation": pd.Series(["=", "="], dtype=pd.StringDtype()),
            "standard_value": pd.Series([1.0, 2.0], dtype="Float64"),
            "standard_units": pd.Series(["NM", "NM"], dtype=pd.StringDtype()),
        }
    )

    activity_data = generate_activity_reports(frame)

    assert isinstance(activity_data, ActivityData)
    assert activity_data.dataset.equals(frame)
    assert activity_data.quality_report.shape[0] == len(frame.columns)
    assert {"row_count", "column_count", "non_null_ratio"} <= activity_data.qc_summary.keys()

