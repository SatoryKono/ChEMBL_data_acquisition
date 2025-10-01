"""Tests for :mod:`library.chembl_assay`."""

from __future__ import annotations

from collections.abc import Iterator
from urllib.parse import parse_qs, urlencode, urlparse

import pandas as pd
from requests import HTTPError, Response

from library.chembl_assay import (
    MAX_TESTITEM_URL_LENGTH,
    TESTITEM_QUERY_FIELDS,
    _split_chunk_for_url,
    get_assays,
    get_testitem,
)
from library.config import ApiCfg


def _response_page(
    items: list[dict[str, str]], next_path: str | None
) -> dict[str, object]:
    """Return a fake paginated response payload."""

    page_meta: dict[str, object] = {"next": next_path}
    return {"assays": items, "page_meta": page_meta}


class DummyClient:
    """Client stub that returns pre-defined responses."""

    def __init__(self, responses: Iterator[dict[str, object]]) -> None:
        self._responses = responses
        self.calls: list[str] = []

    def request_json(
        self, url: str, *, cfg: ApiCfg, timeout: float | None = None
    ) -> dict[str, object]:
        self.calls.append(url)
        try:
            return next(self._responses)
        except StopIteration:  # pragma: no cover - defensive
            raise AssertionError("unexpected extra request") from None


class RecordingTestitemClient:
    """Client stub that records URLs and returns matching molecules."""

    def __init__(self) -> None:
        self.calls: list[str] = []

    def request_json(
        self, url: str, *, cfg: ApiCfg, timeout: float | None = None
    ) -> dict[str, object]:
        self.calls.append(url)
        parsed = urlparse(url)
        params = parse_qs(parsed.query)
        joined = params.get("molecule_chembl_id__in", [""])[0]
        molecules = [
            {"molecule_chembl_id": identifier}
            for identifier in joined.split(",")
            if identifier
        ]
        return {"molecules": molecules, "page_meta": {}}


class LimitedTestitemClient:
    """Client stub that raises HTTP 400 when too many IDs are requested."""

    def __init__(self, max_ids: int) -> None:
        self.max_ids = max_ids
        self.calls: list[str] = []
        self.call_sizes: list[int] = []

    def request_json(
        self, url: str, *, cfg: ApiCfg, timeout: float | None = None
    ) -> dict[str, object]:
        parsed = urlparse(url)
        params = parse_qs(parsed.query)
        joined = params.get("molecule_chembl_id__in", [""])[0]
        identifiers = [identifier for identifier in joined.split(",") if identifier]
        self.calls.append(url)
        self.call_sizes.append(len(identifiers))
        if len(identifiers) > self.max_ids:
            response = Response()
            response.status_code = 400
            response.url = url
            raise HTTPError(
                f"400 Client Error: Bad Request for url: {url}", response=response
            )
        molecules = [
            {"molecule_chembl_id": identifier} for identifier in identifiers
        ]
        return {"molecules": molecules, "page_meta": {}}


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


def test_get_testitem_splits_requests_when_url_would_exceed_limit() -> None:
    """Large batches must be split to keep request URLs under the limit."""

    cfg = ApiCfg()
    base_params: list[tuple[str, str]] = [("format", "json"), ("limit", "1000")]
    if TESTITEM_QUERY_FIELDS:
        base_params.append(("fields", ",".join(TESTITEM_QUERY_FIELDS)))
    base = (
        f"{cfg.chembl_base.rstrip('/')}/molecule.json?{urlencode(base_params)}"
    )
    prefix = f"{base}&molecule_chembl_id__in="
    assert len(prefix) < MAX_TESTITEM_URL_LENGTH

    ids: list[str] = []
    buffer_length = 0
    index = 0
    while True:
        identifier = f"CHEMBL{index:07d}"
        separator = 1 if ids else 0
        candidate_length = buffer_length + separator + len(identifier)
        if len(prefix) + candidate_length > MAX_TESTITEM_URL_LENGTH:
            break
        ids.append(identifier)
        buffer_length = candidate_length
        index += 1

    assert ids, "Base URL must allow at least one identifier"
    ids.extend([f"CHEMBL{index:07d}", f"CHEMBL{index + 1:07d}"])

    splits = list(_split_chunk_for_url(ids, base))
    assert len(splits) >= 2
    for chunk in splits:
        chunk_url = f"{base}&molecule_chembl_id__in={','.join(chunk)}"
        assert len(chunk_url) <= MAX_TESTITEM_URL_LENGTH
    for current, nxt in zip(splits, splits[1:]):
        merged = current + [nxt[0]]
        merged_url = f"{base}&molecule_chembl_id__in={','.join(merged)}"
        assert len(merged_url) > MAX_TESTITEM_URL_LENGTH

    client = RecordingTestitemClient()

    df = get_testitem(
        ids,
        cfg=cfg,
        client=client,
        chunk_size=len(ids),
        page_limit=1000,
    )

    assert len(client.calls) == len(splits)
    assert all(len(url) <= MAX_TESTITEM_URL_LENGTH for url in client.calls)
    assert sorted(df["molecule_chembl_id"].dropna()) == sorted(ids)


def test_get_testitem_splits_chunk_on_http_400() -> None:
    """A 400 response should trigger recursive chunk splitting."""

    cfg = ApiCfg()
    ids = [f"CHEMBL{index:07d}" for index in range(4)]
    client = LimitedTestitemClient(max_ids=2)

    df = get_testitem(
        ids,
        cfg=cfg,
        client=client,
        chunk_size=len(ids),
        page_limit=1000,
    )

    assert sorted(df["molecule_chembl_id"].dropna()) == sorted(ids)
    assert max(client.call_sizes) > 2  # initial oversized request attempted
    assert client.call_sizes.count(2) >= 2  # fallback requests use smaller batches
