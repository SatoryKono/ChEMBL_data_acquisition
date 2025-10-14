"""Integration tests for the standalone activity normalization pipeline."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from library.postprocessing.activity.steps import (
    ACTIVITY_COLUMN_ORDER,
    run_activity_pipeline,
)


class _StubChemblClient:
    """Simple stub returning pre-defined activity pages."""

    def __init__(self, pages: Sequence[dict[str, object]]) -> None:
        self._pages = list(pages)
        self.calls: list[dict[str, object]] = []
        self.timeout = 5.0

    def list_activities(self, **kwargs):
        params = kwargs.get("params", {})
        self.calls.append({"params": dict(params)})
        if not self._pages:
            return {"activities": [], "page_meta": {}}
        return self._pages.pop(0)


def test_run_activity_pipeline__paginates_and_normalizes() -> None:
    pages = [
        {
            "activities": [
                {
                    "activity_id": 20,
                    "assay_chembl_id": "chembl20",
                    "target_chembl_id": "chemblt1",
                    "standard_type": "ic50",
                    "standard_relation": "=",
                    "standard_value": "8.5",
                    "standard_units": "nm",
                }
            ],
            "page_meta": {"next": "token"},
        },
        {
            "activities": [
                {
                    "activity_id": 10,
                    "assay_chembl_id": "chembl10",
                    "target_chembl_id": None,
                    "standard_type": "ec50",
                    "standard_relation": "~=",
                    "standard_value": 2.0,
                    "standard_units": "pm",
                }
            ],
            "page_meta": {},
        },
    ]

    client = _StubChemblClient(pages)
    dataset, correlation, quality = run_activity_pipeline(client=client, page_size=1)

    assert client.calls
    assert list(dataset.columns) == list(ACTIVITY_COLUMN_ORDER)
    assert list(dataset["activity_id"]) == [10, 20]
    assert dataset["standard_units"].tolist() == ["PM", "NM"]
    assert pd.isna(dataset.loc[0, "target_chembl_id"])
    assert dataset.loc[1, "standard_type"] == "IC50"

    assert not quality.empty
    assert "column" in quality.columns
    if not correlation.empty:
        assert list(correlation.index) == list(correlation.columns)

    limited_dataset, _, _ = run_activity_pipeline(limit=1, client=_StubChemblClient(pages=[pages[0]]), page_size=1)
    assert limited_dataset.shape[0] == 1
