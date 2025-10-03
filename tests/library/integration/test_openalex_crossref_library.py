import pytest
import requests

import requests

import pytest

from library.integration import openalex_crossref_library as ocl
from library.common import rate_limiter as rl
from library.config import Config


def test_fetch_openalex_uses_cfg(monkeypatch) -> None:
    """Ensure OpenAlex requests respect configuration parameters."""
    called: dict[str, object] = {}

    def fake_do_request(session, url, delay, timeout=(), **kwargs):
        called["url"] = url
        called["sleep"] = delay
        called["timeout"] = timeout
        return {}, ""

    monkeypatch.setattr("library.clients.openalex._do_request", fake_do_request)

    cfg = Config()
    cfg.api.user_agent = "test@example.com"
    cfg.openalex.base = "https://example.org"
    cfg.openalex.timeout_connect = 1
    cfg.openalex.timeout_read = 2
    cfg.openalex.rps = 2
    cfg.openalex.burst = 5
    cfg.openalex.mailto = "x@y.com"
    limiter = rl.RateLimiter(2)
    ocl.fetch_openalex(requests.Session(), "123", cfg.openalex, limiter)
    assert called["url"] == "https://example.org/works/pmid:123?mailto=x%40y.com"
    assert called["timeout"] == (1, 2)
    assert called["sleep"] == pytest.approx(0.5)


def test_fetch_crossref_uses_cfg(monkeypatch) -> None:
    """Ensure CrossRef requests respect configuration parameters."""
    called: dict[str, object] = {}

    def fake_do_request(session, url, delay, timeout=(), **kwargs):
        called["url"] = url
        called["sleep"] = delay
        called["timeout"] = timeout
        return {}, ""

    monkeypatch.setattr("library.clients.crossref._do_request", fake_do_request)

    cfg = Config()
    cfg.api.user_agent = "test@example.com"
    cfg.crossref.base = "https://cr.example.org"
    cfg.crossref.timeout_connect = 1
    cfg.crossref.timeout_read = 2
    cfg.crossref.rps = 4
    cfg.crossref.burst = 5
    cfg.crossref.mailto = "z@e.com"
    limiter = rl.RateLimiter(4)
    ocl.fetch_crossref(requests.Session(), "10.1/abc", cfg.crossref, limiter)
    assert called["url"] == "https://cr.example.org/works/10.1%2Fabc?mailto=z%40e.com"
    assert called["timeout"] == (1, 2)
    assert called["sleep"] == pytest.approx(0.25)


def test_rate_limiter_shared(monkeypatch) -> None:
    """Limiter should throttle subsequent calls when reused."""
    delays: list[float] = []

    def fake_sleep(delay: float) -> None:
        delays.append(delay)

    monkeypatch.setattr(rl, "sleep", fake_sleep)
    monkeypatch.setattr(
        "library.clients.openalex._do_request", lambda *a, **k: ({}, "")
    )

    cfg = Config()
    cfg.api.user_agent = "test@example.com"
    cfg.openalex.rps = 1
    cfg.openalex.mailto = "x@y.com"
    limiter = rl.RateLimiter(1)
    session = requests.Session()
    ocl.fetch_openalex(session, "1", cfg.openalex, limiter)
    ocl.fetch_openalex(session, "2", cfg.openalex, limiter)
    assert delays and delays[0] == pytest.approx(1.0, rel=0.1)
