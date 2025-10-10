"""Unit tests for :mod:`library.clients.chembl`."""

from __future__ import annotations

import json
import socket
from dataclasses import dataclass
from http.client import HTTPConnection
from urllib.parse import urlsplit

import pytest
import requests
from urllib3.connectionpool import HTTPSConnectionPool
from urllib3.exceptions import MaxRetryError, NameResolutionError, ReadTimeoutError

from library.clients import ChemblClient
from library.clients.chembl import _backoff_delay, _normalise_request_exception
from library.config import ApiCfg, RetryCfg


@dataclass
class _StubResponse:
    """Minimal stand-in for :class:`requests.Response`."""

    url: str
    status_code: int
    payload: dict[str, object] | None = None
    headers: dict[str, str] | None = None

    def __post_init__(self) -> None:
        self.headers = self.headers or {}

    def __enter__(self) -> _StubResponse:  # pragma: no cover - context API
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


class _TimeoutSession:
    """Raise a timeout on the primary URL before succeeding on fallback."""

    def __init__(self, primary: str, fallback: str, payload: dict[str, object]) -> None:
        self.primary_url = primary
        self.fallback_url = fallback
        self.payload = payload
        self.calls: list[str] = []

    def get(self, url: str, timeout: object) -> _StubResponse:
        del timeout
        self.calls.append(url)
        if url == self.primary_url:
            raise requests.ReadTimeout(f"timeout for {url}")
        if url == self.fallback_url:
            return _StubResponse(url, 200, self.payload)
        raise AssertionError(f"Unexpected URL requested: {url}")

    def close(self) -> None:  # pragma: no cover - compatibility shim
        return None


class _MaxRetryTimeoutSession:
    """Raise a wrapped timeout resembling urllib3's ``MaxRetryError``."""

    def __init__(self, primary: str, fallback: str) -> None:
        self.primary_url = primary
        self.fallback_url = fallback
        self.calls: list[str] = []

    def get(self, url: str, timeout: object) -> _StubResponse:
        del timeout
        self.calls.append(url)
        host = urlsplit(url).netloc or "example.test"
        pool = HTTPSConnectionPool(host)
        reason = ReadTimeoutError(pool, url, "Read timed out.")
        raise requests.ConnectionError(MaxRetryError(pool, url, reason))

    def close(self) -> None:  # pragma: no cover - compatibility shim
        return None


class _Always404Session:
    """Return 404 for both primary and fallback URLs."""

    def __init__(self, primary: str, fallback: str) -> None:
        self.primary_url = primary
        self.fallback_url = fallback
        self.calls: list[str] = []

    def get(self, url: str, timeout: object) -> _StubResponse:
        del timeout
        self.calls.append(url)
        return _StubResponse(url, 404)

    def close(self) -> None:  # pragma: no cover - compatibility shim
        return None


class _EmptyResponseSession:
    """Return an empty payload that triggers JSON parsing errors."""

    def __init__(self, primary: str) -> None:
        self.primary_url = primary
        self.calls: list[str] = []

    def get(self, url: str, timeout: object) -> _StubResponse:
        del timeout
        self.calls.append(url)
        return _StubResponse(url, 200, payload=None)

    def close(self) -> None:  # pragma: no cover - compatibility shim
        return None


class _InvalidJSONSession:
    """Return a response object raising ``JSONDecodeError``."""

    def __init__(self, primary: str) -> None:
        self.primary_url = primary
        self.calls: list[str] = []

    def get(self, url: str, timeout: object) -> _StubResponse:
        del timeout
        self.calls.append(url)

        class _InvalidResponse(_StubResponse):
            def json(self) -> dict[str, object]:  # type: ignore[override]
                raise json.JSONDecodeError("Expecting value", "", 0)

        return _InvalidResponse(url, 200, payload={})

    def close(self) -> None:  # pragma: no cover - compatibility shim
        return None


class _NameResolutionSession:
    """Raise ``NameResolutionError`` to emulate DNS failures."""

    def __init__(self, primary: str) -> None:
        self.primary_url = primary
        self.calls: list[str] = []

    def get(self, url: str, timeout: object) -> _StubResponse:
        del timeout
        self.calls.append(url)
        host = urlsplit(url).hostname or "example.test"
        connection = HTTPConnection(host)
        reason = socket.gaierror(-2, "Name or service not known")
        raise requests.ConnectionError(NameResolutionError(host, connection, reason))

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


@pytest.mark.unit
def test_request_json__falls_back_after_timeout() -> None:
    """Connection timeouts trigger the extensionless fallback."""

    base = "https://example.test/chembl/api/data"
    query = "format=json&assay_chembl_id__in=CHEMBL1&limit=1"
    primary_url = f"{base}/assay.json?{query}"
    fallback_url = f"{base}/assay?{query}"
    payload = {"assays": [{"assay_chembl_id": "CHEMBL1"}]}
    session = _TimeoutSession(primary_url, fallback_url, payload)
    client = ChemblClient(session=session)
    cfg = ApiCfg(chembl_base=base, timeout_read=5.0)

    result = client.request_json(primary_url, cfg=cfg)

    assert result == payload
    assert session.calls == [primary_url, fallback_url]


@pytest.mark.unit
def test_request_json__raises_read_timeout_for_max_retry_chain() -> None:
    """Connection errors caused by read timeouts surface as ``ReadTimeout``."""

    base = "https://example.test/chembl/api/data"
    query = "format=json&assay_chembl_id__in=CHEMBL1&limit=1"
    primary_url = f"{base}/assay.json?{query}"
    fallback_url = f"{base}/assay?{query}"
    session = _MaxRetryTimeoutSession(primary_url, fallback_url)
    client = ChemblClient(session=session)
    cfg = ApiCfg(chembl_base=base, timeout_read=5.0, retries=0)

    with pytest.raises(requests.ReadTimeout) as excinfo:
        client.request_json(primary_url, cfg=cfg)

    assert "Read timed out" in str(excinfo.value)
    assert session.calls == [primary_url, fallback_url]


@pytest.mark.unit
def test_request_json__raises_for_persistent_404() -> None:
    """The client surfaces HTTP errors when fallback also returns 404."""

    base = "https://example.test/chembl/api/data"
    query = "format=json&assay_chembl_id__in=CHEMBL1&limit=1"
    primary_url = f"{base}/assay.json?{query}"
    fallback_url = f"{base}/assay?{query}"
    session = _Always404Session(primary_url, fallback_url)
    client = ChemblClient(session=session)
    cfg = ApiCfg(chembl_base=base, timeout_read=5.0, retries=0)

    with pytest.raises(requests.HTTPError):
        client.request_json(primary_url, cfg=cfg)

    assert session.calls == [primary_url, fallback_url]


@pytest.mark.unit
def test_request_json__raises_for_empty_response() -> None:
    """Empty payloads raise ``ValueError`` signalling invalid JSON."""

    base = "https://example.test/chembl/api/data"
    primary_url = f"{base}/assay.json?format=json"
    session = _EmptyResponseSession(primary_url)
    client = ChemblClient(session=session)
    cfg = ApiCfg(chembl_base=base, timeout_read=5.0, retries=0)

    with pytest.raises(ValueError):
        client.request_json(primary_url, cfg=cfg)

    assert session.calls == [primary_url]


@pytest.mark.unit
def test_request_json__raises_for_invalid_json() -> None:
    """Malformed JSON surfaces as a ``ValueError`` to the caller."""

    base = "https://example.test/chembl/api/data"
    primary_url = f"{base}/assay.json?format=json"
    session = _InvalidJSONSession(primary_url)
    client = ChemblClient(session=session)
    cfg = ApiCfg(chembl_base=base, timeout_read=5.0, retries=0)

    with pytest.raises(ValueError):
        client.request_json(primary_url, cfg=cfg)

    assert session.calls == [primary_url]


@pytest.mark.unit
def test_request_json__fails_fast_on_name_resolution_error() -> None:
    """DNS failures are not retried and raise ``ConnectionError``."""

    base = "https://example.test/chembl/api/data"
    primary_url = f"{base}/assay.json?format=json"
    session = _NameResolutionSession(primary_url)
    client = ChemblClient(session=session)
    cfg = ApiCfg(chembl_base=base, timeout_read=5.0, retries=3)

    with pytest.raises(requests.ConnectionError):
        client.request_json(primary_url, cfg=cfg)

    assert session.calls == [primary_url]


@pytest.mark.unit
def test_backoff_delay__deterministic_jitter() -> None:
    api_cfg = ApiCfg(backoff_factor=0.5)
    jitter_one = RetryCfg(jitter_seed=11).build_jitter()
    jitter_two = RetryCfg(jitter_seed=11).build_jitter()

    assert jitter_one is not None
    assert jitter_two is not None

    delays_one = [
        _backoff_delay(attempt, api_cfg, None, jitter=jitter_one)
        for attempt in range(1, 4)
    ]
    delays_two = [
        _backoff_delay(attempt, api_cfg, None, jitter=jitter_two)
        for attempt in range(1, 4)
    ]

    assert delays_one == delays_two


@pytest.mark.unit
def test_normalise_request_exception__read_timeout_message() -> None:
    connection_error = requests.ConnectionError(
        "HTTPSConnectionPool(host='www.ebi.ac.uk', port=443): Read timed out."
    )

    normalised = _normalise_request_exception(connection_error)

    assert normalised is not connection_error
    assert isinstance(normalised, requests.ReadTimeout)
    assert "read" in str(normalised).lower()


@pytest.mark.unit
def test_normalise_request_exception__preserves_request_metadata() -> None:
    timeout = requests.ReadTimeout("Read timed out.")
    timeout.request = object()  # type: ignore[attr-defined]
    timeout.response = object()  # type: ignore[attr-defined]
    exc = requests.ConnectionError(timeout)

    normalised = _normalise_request_exception(exc)

    assert isinstance(normalised, requests.ReadTimeout)
    assert normalised.request is timeout.request
    assert normalised.response is timeout.response
