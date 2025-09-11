import pytest
import requests

import library.rate_limiter as rl
from library import openalex_crossref_library as ocl
from library.config import Config, CrossRefCfg, OpenAlexCfg


def test_fetch_openalex_uses_cfg(monkeypatch) -> None:
    """Ensure OpenAlex requests respect configuration parameters."""
    called: dict[str, object] = {}

    def fake_do_request(session, url, delay, timeout=(), **kwargs):
        called["url"] = url
        called["sleep"] = delay
        called["timeout"] = timeout
        return {}, ""

    monkeypatch.setattr("library.pubmed_library._do_request", fake_do_request)

    cfg = Config(
        openalex=OpenAlexCfg(
            base="https://example.org",
            timeout_connect=1,
            timeout_read=2,
            rps=2,
            burst=5,
            mailto="x@y.com",
        )
    )
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

    monkeypatch.setattr("library.pubmed_library._do_request", fake_do_request)

    cfg = Config(
        crossref=CrossRefCfg(
            base="https://cr.example.org",
            timeout_connect=1,
            timeout_read=2,
            rps=4,
            burst=5,
            mailto="z@e.com",
        )
    )
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
    monkeypatch.setattr("library.pubmed_library._do_request", lambda *a, **k: ({}, ""))

    cfg = Config(openalex=OpenAlexCfg(rps=1, mailto="x@y.com"))
    limiter = rl.RateLimiter(1)
    session = requests.Session()
    ocl.fetch_openalex(session, "1", cfg.openalex, limiter)
    ocl.fetch_openalex(session, "2", cfg.openalex, limiter)
    assert delays and delays[0] == pytest.approx(1.0, rel=0.1)
