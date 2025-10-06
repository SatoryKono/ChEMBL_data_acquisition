"""Unit tests for :mod:`library.pipelines.assay.chembl_assay`."""

from __future__ import annotations

from typing import Any

import pytest
import requests
from requests import HTTPError, Response

from library.config import ApiCfg
from library.pipelines.assay.chembl_assay import get_assays


class _StubClient:
    """Simple stub mimicking the :class:`ChemblClient` interface."""

    def __init__(self, responders: list[dict[str, Any] | str]) -> None:
        self._responders = responders
        self.calls: list[str] = []

    def request_json(self, url: str, *, cfg: ApiCfg, timeout: float | None) -> dict[str, Any]:
        self.calls.append(url)
        next_response = self._responders.pop(0)
        if isinstance(next_response, str):
            if next_response == "HTTP404":
                response = Response()
                response.status_code = 404
                response.url = url
                raise HTTPError(response=response)
            if next_response == "TIMEOUT":
                raise requests.Timeout("simulated timeout")
        if isinstance(next_response, dict):
            return next_response
        raise AssertionError(f"Unexpected responder type: {type(next_response)!r}")


@pytest.mark.unit
def test_get_assays__splits_chunk_on_404() -> None:
    """Chunk-level 404 errors fall back to per-ID requests."""

    responders = [
        "HTTP404",
        {"assays": [{"assay_chembl_id": "CHEMBL1"}], "page_meta": {}},
        {"assays": [], "page_meta": {}},
    ]
    client = _StubClient(responders)
    cfg = ApiCfg()

    df = get_assays(["CHEMBL1", "CHEMBL2"], cfg=cfg, client=client, chunk_size=2)

    assert list(df["assay_chembl_id"]) == ["CHEMBL1"]
    assert any("assay_chembl_id__in" in call for call in client.calls)
    assert any("assay_chembl_id=CHEMBL1" in call for call in client.calls)
    assert any("assay_chembl_id=CHEMBL2" in call for call in client.calls)


@pytest.mark.unit
def test_get_assays__detail_fallback_on_single_404() -> None:
    """Single-item 404 fallbacks query the detail endpoint."""

    responders = [
        "HTTP404",
        "HTTP404",
        {"assay": {"assay_chembl_id": "CHEMBL1"}},
        "HTTP404",
        {"assays": []},
    ]
    client = _StubClient(responders)
    cfg = ApiCfg()

    df = get_assays(["CHEMBL1", "CHEMBL2"], cfg=cfg, client=client, chunk_size=2)

    assert list(df["assay_chembl_id"]) == ["CHEMBL1"]
    assert any("/assay/CHEMBL1" in call for call in client.calls)


@pytest.mark.unit
def test_get_assays__splits_chunk_on_timeout() -> None:
    """Timeouts trigger chunk splitting retries."""

    responders = [
        "TIMEOUT",
        {"assays": [{"assay_chembl_id": "CHEMBL1"}], "page_meta": {}},
        {"assays": [{"assay_chembl_id": "CHEMBL2"}], "page_meta": {}},
    ]
    client = _StubClient(responders)
    cfg = ApiCfg()

    df = get_assays(["CHEMBL1", "CHEMBL2"], cfg=cfg, client=client, chunk_size=2)

    assert sorted(df["assay_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
    assert any("limit=2" in call for call in client.calls)
    assert sum(1 for call in client.calls if "limit=1" in call) == 2


@pytest.mark.unit
def test_get_assays__preserves_data_segment_in_pagination() -> None:
    """Relative ``page_meta.next`` links keep the ``/data`` segment."""

    responders = [
        {
            "assays": [{"assay_chembl_id": "CHEMBL1"}],
            "page_meta": {"next": "assay.json?offset=1&limit=1"},
        },
        {"assays": [{"assay_chembl_id": "CHEMBL2"}], "page_meta": {}},
    ]
    client = _StubClient(responders)
    cfg = ApiCfg()

    df = get_assays(["CHEMBL1"], cfg=cfg, client=client, chunk_size=1)

    assert sorted(df["assay_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
    assert any(
        call.startswith("https://www.ebi.ac.uk/chembl/api/data/assay.json?offset=1")
        for call in client.calls[1:]
    )
