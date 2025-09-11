import requests
import pytest

from library import openalex_crossref_library as ocl
from library.config import OpenAlexCfg, CrossRefCfg


def test_fetch_openalex_uses_cfg(monkeypatch) -> None:
    """Ensure OpenAlex requests respect configuration parameters."""
    called: dict[str, object] = {}

    def fake_do_request(session, url, delay, timeout=(), **kwargs):
        called["url"] = url
        called["sleep"] = delay
        called["timeout"] = timeout
        return {}, ""

    monkeypatch.setattr("library.pubmed_library._do_request", fake_do_request)
    rps: dict[str, float] = {}

    def fake_get_limiter(name: str, rps_val: float, burst: int | None = None):
        rps["value"] = rps_val

        class Dummy:
            def acquire(self) -> None:
                pass

        return Dummy()

    monkeypatch.setattr("library.pubmed_library.get_limiter", fake_get_limiter)

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
    assert rps["value"] == pytest.approx(2)


def test_fetch_crossref_uses_cfg(monkeypatch) -> None:
    """Ensure CrossRef requests respect configuration parameters."""
    called: dict[str, object] = {}

    def fake_do_request(session, url, delay, timeout=(), **kwargs):
        called["url"] = url
        called["sleep"] = delay
        called["timeout"] = timeout
        return {}, ""

    monkeypatch.setattr("library.pubmed_library._do_request", fake_do_request)
    rps: dict[str, float] = {}

    def fake_get_limiter(name: str, rps_val: float, burst: int | None = None):
        rps["value"] = rps_val

        class Dummy:
            def acquire(self) -> None:
                pass

        return Dummy()

    monkeypatch.setattr("library.pubmed_library.get_limiter", fake_get_limiter)

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
    assert rps["value"] == pytest.approx(4)
