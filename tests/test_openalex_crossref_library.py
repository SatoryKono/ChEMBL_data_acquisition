import requests
import pytest

from library import openalex_crossref_library as ocl
from library.config import OpenAlexCfg, CrossRefCfg


def test_fetch_openalex_uses_cfg(monkeypatch) -> None:
    """Ensure OpenAlex requests respect configuration parameters."""
    called: dict[str, object] = {}

    def fake_do_request(session, url, sleep, timeout=(), **kwargs):
        called["url"] = url
        called["sleep"] = sleep
        called["timeout"] = timeout
        return {}, ""

    monkeypatch.setattr("library.pubmed_library._do_request", fake_do_request)
    sleeps: list[float] = []
    monkeypatch.setattr("library.pubmed_library.time.sleep", lambda s: sleeps.append(s))

    cfg = OpenAlexCfg(
        base="https://example.org",
        timeout_connect=1,
        timeout_read=2,
        rps=2,
        burst=5,
        mailto="x@y.com",
    )
    ocl.fetch_openalex(requests.Session(), "123", cfg)
    assert called["url"] == "https://example.org/works/pmid:123?mailto=x%40y.com"
    assert called["timeout"] == (1, 2)
    assert sleeps and sleeps[0] == pytest.approx(0.5)


def test_fetch_crossref_uses_cfg(monkeypatch) -> None:
    """Ensure CrossRef requests respect configuration parameters."""
    called: dict[str, object] = {}

    def fake_do_request(session, url, sleep, timeout=(), **kwargs):
        called["url"] = url
        called["sleep"] = sleep
        called["timeout"] = timeout
        return {}, ""

    monkeypatch.setattr("library.pubmed_library._do_request", fake_do_request)
    sleeps: list[float] = []
    monkeypatch.setattr("library.pubmed_library.time.sleep", lambda s: sleeps.append(s))

    cfg = CrossRefCfg(
        base="https://cr.example.org",
        timeout_connect=1,
        timeout_read=2,
        rps=4,
        burst=5,
        mailto="z@e.com",
    )
    ocl.fetch_crossref(requests.Session(), "10.1/abc", cfg)
    assert called["url"] == "https://cr.example.org/works/10.1%2Fabc?mailto=z%40e.com"
    assert called["timeout"] == (1, 2)
    assert sleeps and sleeps[0] == pytest.approx(0.25)
