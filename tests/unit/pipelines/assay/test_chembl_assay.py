"""Unit tests for :mod:`library.pipelines.assay.chembl_assay`."""

from __future__ import annotations

from typing import Any, Callable

import pytest
from requests import ConnectionError, HTTPError, ReadTimeout, Response
from urllib.parse import parse_qs, urlparse

from library.config import ApiCfg
from library.pipelines.assay.chembl_assay import get_assays, get_testitem


class _StubClient:
    """Simple stub mimicking the :class:`ChemblClient` interface."""

    def __init__(self, responders: list[dict[str, Any] | str]) -> None:
        self._responders = responders
        self.calls: list[str] = []

    def request_json(self, url: str, *, cfg: ApiCfg, timeout: float | None) -> dict[str, Any]:
        self.calls.append(url)
        next_response = self._responders.pop(0)
        if isinstance(next_response, str) and next_response == "HTTP404":
            response = Response()
            response.status_code = 404
            response.url = url
            raise HTTPError(response=response)
        if isinstance(next_response, dict):
            return next_response
        if isinstance(next_response, Exception):
            raise next_response
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


@pytest.mark.unit
def test_get_assays__recovers_from_request_exception() -> None:
    """Network failures fall back to per-ID requests."""

    responders = [
        ReadTimeout("timeout"),
        {"assays": [{"assay_chembl_id": "CHEMBL1"}], "page_meta": {}},
        {"assays": [{"assay_chembl_id": "CHEMBL2"}], "page_meta": {}},
    ]
    client = _StubClient(responders)
    cfg = ApiCfg()

    df = get_assays(["CHEMBL1", "CHEMBL2"], cfg=cfg, client=client, chunk_size=2)

    assert sorted(df["assay_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
    assert any("assay_chembl_id__in" in call for call in client.calls)
    assert sum("assay_chembl_id=CHEMBL1" in call for call in client.calls) == 1
    assert sum("assay_chembl_id=CHEMBL2" in call for call in client.calls) == 1


@pytest.mark.unit
def test_get_testitem__splits_chunk_on_timeout() -> None:
    """Timeout failures trigger chunk splitting and retries."""

    def _timeout_exc() -> ConnectionError:
        return ConnectionError(
            "HTTPSConnectionPool(host='www.ebi.ac.uk', port=443): Read timed out."
        )

    payloads: dict[tuple[str, ...], dict[str, Any] | Callable[[], Exception]] = {
        ("CHEMBL1", "CHEMBL2", "CHEMBL3"): _timeout_exc,
        ("CHEMBL1",): {
            "molecules": [
                {"molecule_chembl_id": "CHEMBL1", "pref_name": "One"},
            ],
            "page_meta": {},
        },
        ("CHEMBL2", "CHEMBL3"): {
            "molecules": [
                {"molecule_chembl_id": "CHEMBL2", "pref_name": "Two"},
                {"molecule_chembl_id": "CHEMBL3", "pref_name": "Three"},
            ],
            "page_meta": {},
        },
    }

    class _TimeoutStub:
        def __init__(
            self,
            mapping: dict[tuple[str, ...], dict[str, Any] | Callable[[], Exception]],
        ) -> None:
            self._mapping = mapping
            self.calls: list[str] = []

        def request_json(self, url: str, *, cfg: ApiCfg, timeout: float | None) -> dict[str, Any]:
            del cfg, timeout
            self.calls.append(url)
            parsed = urlparse(url)
            ids_param = parse_qs(parsed.query).get("molecule_chembl_id__in", [""])[0]
            identifiers = tuple(filter(None, ids_param.split(",")))
            response = self._mapping[identifiers]
            if callable(response):
                response = response()
            if isinstance(response, Exception):
                raise response
            return response

    client = _TimeoutStub(payloads)
    cfg = ApiCfg()

    df = get_testitem(
        ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
        cfg=cfg,
        client=client,  # type: ignore[arg-type]
        chunk_size=3,
    )

    assert sorted(df["molecule_chembl_id"].dropna().tolist()) == [
        "CHEMBL1",
        "CHEMBL2",
        "CHEMBL3",
    ]
    assert any(
        "molecule_chembl_id__in=CHEMBL1,CHEMBL2,CHEMBL3" in call for call in client.calls
    )
    assert any("molecule_chembl_id__in=CHEMBL1" in call for call in client.calls)
    assert any(
        "molecule_chembl_id__in=CHEMBL2,CHEMBL3" in call for call in client.calls
    )
def test_get_assays__single_timeout_falls_back_to_detail_endpoint() -> None:
    """Single-item timeouts fall back to the detail endpoint."""

    responders = [
        ReadTimeout("timeout"),
        ReadTimeout("timeout"),
        {"assay": {"assay_chembl_id": "CHEMBL1", "sequence": "AA"}},
        {"assays": [{"assay_chembl_id": "CHEMBL2", "sequence": "BB"}], "page_meta": {}},
    ]
    client = _StubClient(responders)
    cfg = ApiCfg()

    df = get_assays(
        ["CHEMBL1", "CHEMBL2"],
        cfg=cfg,
        client=client,
        chunk_size=2,
        require_variant_sequence=True,
    )

    assert sorted(df["assay_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
    assert any("assay_chembl_id__in" in call for call in client.calls)
    assert any("/assay/CHEMBL1" in call for call in client.calls)
