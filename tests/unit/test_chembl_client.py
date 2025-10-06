"""Unit tests for :mod:`library.clients.chembl`."""

from __future__ import annotations

from dataclasses import dataclass

import pytest
import requests

from library.clients import ChemblClient
from library.config import ApiCfg


@dataclass
class _StubResponse:
    """Minimal stand-in for :class:`requests.Response`."""

    url: str
    status_code: int
    payload: dict[str, object] | None = None
    headers: dict[str, str] | None = None

    def __post_init__(self) -> None:
        self.headers = self.headers or {}

    def __enter__(self) -> "_StubResponse":  # pragma: no cover - context API
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb,
    ) -> bool:  # pragma: no cover - context API
        return False

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(
                f"{self.status_code} error for {self.url}", response=self
            )

    def json(self) -> dict[str, object]:
        if self.payload is None:
            raise ValueError("no JSON payload available")
        return self.payload


class _StubSession:
    """Return deterministic responses for ChemblClient tests."""

    def __init__(self, primary: str, fallback: str, payload: dict[str, object]) -> None:
        self.primary_url = primary
        self.fallback_url = fallback
        self.payload = payload
        self.calls: list[str] = []

    def get(self, url: str, timeout: object) -> _StubResponse:
        self.calls.append(url)
        if url == self.primary_url:
            return _StubResponse(url, 404)
        if url == self.fallback_url:
            return _StubResponse(url, 200, self.payload)
        raise AssertionError(f"Unexpected URL requested: {url}")

    def close(self) -> None:  # pragma: no cover - compatibility shim
        return None


@pytest.mark.unit
def test_request_json__falls_back_to_extensionless_endpoint() -> None:
    """The client retries with an extensionless URL when a 404 is returned."""

    base = "https://example.test/chembl/api/data"
    query = "format=json&assay_chembl_id__in=CHEMBL1&limit=1"
    primary_url = f"{base}/assay.json?{query}"
    fallback_url = f"{base}/assay?{query}"
    payload = {"assays": [{"assay_chembl_id": "CHEMBL1"}]}
    session = _StubSession(primary_url, fallback_url, payload)
    client = ChemblClient(session=session)
    cfg = ApiCfg(chembl_base=base, timeout_read=5.0)

    result = client.request_json(primary_url, cfg=cfg)

    assert result == payload
    assert session.calls == [primary_url, fallback_url]

    cached = client.request_json(primary_url, cfg=cfg)
    assert cached == payload
    assert session.calls == [primary_url, fallback_url]
