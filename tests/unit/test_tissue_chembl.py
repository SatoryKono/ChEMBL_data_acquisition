"""Unit tests for :func:`library.pipelines.tissue.chembl.get_tissues`."""

from __future__ import annotations

from collections import deque
from typing import Any, Iterable

import pandas as pd
import pytest

from library.config import ApiCfg
from library.pipelines.tissue.chembl import TISSUE_BASE_COLUMNS, get_tissues


class MockChemblClient:
    """Minimal stand-in for :class:`~library.clients.ChemblClient` used in tests."""

    def __init__(self, payloads: Iterable[dict[str, Any]]):
        self._payloads: deque[dict[str, Any]] = deque(payloads)
        self.calls: list[tuple[str, float | None]] = []

    def request_json(
        self,
        url: str,
        *,
        cfg: ApiCfg,
        timeout: float | None = None,
    ) -> dict[str, Any]:
        self.calls.append((url, timeout))
        if not self._payloads:
            raise AssertionError("No payload queued for request")
        return self._payloads.popleft()


@pytest.mark.unit
def test_get_tissues__chunked_requests() -> None:
    cfg = ApiCfg(chembl_base="https://example.test/custom", timeout_read=12.5)
    client = MockChemblClient(
        [
            {
                "tissues": [
                    {
                        "tissue_chembl_id": "TIS_1",
                        "pref_name": "Alpha",
                        "uberon_id": "UBERON:1",
                    }
                ],
                "page_meta": {
                    "next": "https://example.test/custom/tissue.json?format=json&offset=1"
                },
            },
            {
                "tissues": [
                    {
                        "tissue_chembl_id": "TIS_2",
                        "pref_name": "Beta",
                        "uberon_id": "UBERON:2",
                    }
                ],
                "page_meta": {},
            },
            {
                "tissues": [
                    {
                        "tissue_chembl_id": "TIS_3",
                        "pref_name": "Gamma",
                        "uberon_id": "UBERON:3",
                    },
                    {
                        "tissue_chembl_id": "TIS_4",
                        "pref_name": "Delta",
                        "uberon_id": "UBERON:4",
                    },
                ],
                "page_meta": {},
            },
        ]
    )

    df = get_tissues(
        ["TIS_1", "TIS_2", "TIS_3", "TIS_4"],
        cfg=cfg,
        client=client,
        chunk_size=2,
    )

    expected_urls = [
        "https://example.test/custom/tissue.json?format=json"
        "&tissue_chembl_id__in=TIS_1,TIS_2&limit=2",
        "https://example.test/custom/tissue.json?format=json&offset=1",
        "https://example.test/custom/tissue.json?format=json"
        "&tissue_chembl_id__in=TIS_3,TIS_4&limit=2",
    ]
    assert client.calls == [(url, 12.5) for url in expected_urls]
    assert df["tissue_chembl_id"].tolist() == ["TIS_1", "TIS_2", "TIS_3", "TIS_4"]


@pytest.mark.unit
def test_get_tissues__normalises_missing_columns() -> None:
    cfg = ApiCfg(chembl_base="https://example.test/base", timeout_read=5.0)
    client = MockChemblClient(
        [
            {
                "tissues": [
                    {
                        "tissue_chembl_id": "TIS_MISSING",
                        "pref_name": None,
                        "uberon_id": "",
                        "caloha_id": None,
                    }
                ],
                "page_meta": {},
            }
        ]
    )

    df = get_tissues(
        ["TIS_MISSING"],
        cfg=cfg,
        client=client,
        chunk_size=5,
    )

    assert df.columns.tolist() == TISSUE_BASE_COLUMNS
    row = df.iloc[0]
    assert row["tissue_chembl_id"] == "TIS_MISSING"
    assert pd.isna(row["pref_name"])
    assert row["uberon_id"] == ""
    assert pd.isna(row["efo_id"])
    assert pd.isna(row["bto_id"])
    assert pd.isna(row["caloha_id"])


@pytest.mark.unit
def test_get_tissues__skips_invalid_ids() -> None:
    cfg = ApiCfg(chembl_base="https://example.test/base", timeout_read=5.0)
    client = MockChemblClient([])

    df = get_tissues(
        ["", "#N/A"],
        cfg=cfg,
        client=client,
        chunk_size=3,
    )

    assert client.calls == []
    assert df.empty
    assert df.columns.tolist() == TISSUE_BASE_COLUMNS

