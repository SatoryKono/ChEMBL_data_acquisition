"""Tests for :mod:`library.chembl_assay`."""

from __future__ import annotations

from typing import Iterator

import pandas as pd

from library.chembl_assay import get_assays
from library.config import ApiCfg


def _response_page(items: list[dict[str, str]], next_path: str | None) -> dict[str, object]:
    """Return a fake paginated response payload."""

    page_meta: dict[str, object] = {"next": next_path}
    return {"assays": items, "page_meta": page_meta}


class DummyClient:
    """Client stub that returns pre-defined responses."""

    def __init__(self, responses: Iterator[dict[str, object]]) -> None:
        self._responses = responses
        self.calls: list[str] = []

    def request_json(self, url: str, *, cfg: ApiCfg, timeout: float | None = None) -> dict[str, object]:
        self.calls.append(url)
        try:
            return next(self._responses)
        except StopIteration:  # pragma: no cover - defensive
            raise AssertionError("unexpected extra request") from None


def test_get_assays_fetches_all_pages() -> None:
    """Pagination should be followed when ChEMBL truncates responses."""

    cfg = ApiCfg()
    ids = ["A1", "A2", "A3", "A4"]
    next_path = "/chembl/api/data/assay.json?format=json&assay_chembl_id__in=A1,A2,A3,A4&limit=4&offset=2"
    responses = iter(
        [
            _response_page(
                [
                    {"assay_chembl_id": "A1", "document_chembl_id": "D1"},
                    {"assay_chembl_id": "A2", "document_chembl_id": "D2"},
                ],
                next_path,
            ),
            _response_page(
                [
                    {"assay_chembl_id": "A3", "document_chembl_id": "D3"},
                    {"assay_chembl_id": "A4", "document_chembl_id": "D4"},
                ],
                None,
            ),
        ]
    )
    client = DummyClient(responses)

    df = get_assays(ids, cfg=cfg, client=client, chunk_size=4)

    assert len(client.calls) == 2
    assert "limit=4" in client.calls[0]
    assert "offset=2" in client.calls[1]
    assert isinstance(df, pd.DataFrame)
    assert list(df["assay_chembl_id"]) == ids
