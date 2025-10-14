"""Unit tests covering timestamp derivation in the assay pipeline."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

from library.postprocessing.assay.steps import AssayFetcher, fetch_normalize_assay


def test_timestamp_and_year_are_derived(monkeypatch):
    records = [
        {
            "assay_chembl_id": "CHEMBL1",
            "assay_type": "BINDING",
            "assay_test_type": "IN VITRO",
            "target_chembl_id": "CHEMBL10",
            "created_on": "2020-05-01T12:00:00",
            "updated_on": None,
            "assay_strain": None,
            "assay_group": None,
            "accession": None,
        },
        {
            "assay_chembl_id": "CHEMBL2",
            "assay_type": "BINDING",
            "assay_test_type": "IN VITRO",
            "target_chembl_id": "CHEMBL20",
            "created_on": None,
            "updated_on": "2021-03-15T08:30:00",
            "assay_strain": None,
            "assay_group": None,
            "accession": None,
        },
    ]

    def fake_fetch(self, limit=None):  # type: ignore[unused-argument]
        return list(records)

    monkeypatch.setattr(AssayFetcher, "fetch", fake_fetch, raising=False)

    frame = fetch_normalize_assay(dictionary_dir=Path("non-existent"))

    assert list(frame["assay_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
    assert frame["timestamp_utc"].dtype == "datetime64[ns, UTC]"

    expected = [
        datetime(2020, 5, 1, 12, 0, tzinfo=timezone.utc),
        datetime(2021, 3, 15, 8, 30, tzinfo=timezone.utc),
    ]
    assert frame["timestamp_utc"].tolist() == expected
    assert frame["year"].tolist() == [2020, 2021]
    assert str(frame["year"].dtype) == "Int64"
