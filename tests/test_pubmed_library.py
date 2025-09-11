"""Tests for :mod:`library.pubmed_library`."""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import pytest
import requests

import library.rate_limiter as rl
from library import pubmed_library as pl
from library.config import (
    Config,
    CrossRefCfg,
    OpenAlexCfg,
    PubMedCfg,
    SemanticScholarCfg,
)

DATA_DIR = Path(__file__).parent / "data"


def test_read_pmids() -> None:
    path = DATA_DIR / "pmids.csv"
    pmids = pl.read_pmids(path)
    assert pmids == ["1", "2"]


def test_read_pmids_missing_column(tmp_path: Path) -> None:
    """Ensure ``read_pmids`` validates presence of the PMID column."""
    path = tmp_path / "bad.csv"
    path.write_text("ID\n1\n")
    with pytest.raises(ValueError, match="PMID"):
        pl.read_pmids(path)


class DummyResponse:
    def __init__(
        self, status: int, text: str = "", json_data: dict[str, Any] | None = None
    ) -> None:
        self.status_code = status
        self._text = text
        self._json = json_data

    def __enter__(self) -> DummyResponse:
        return self

    def __exit__(
        self, exc_type, exc, tb
    ) -> None:  # pragma: no cover - no cleanup needed
        return None

    @property
    def text(self) -> str:
        return self._text

    def json(self) -> dict[str, Any]:
        if self._json is None:
            raise ValueError("no json")
        return self._json


class DummySession:
    def __init__(self, response: DummyResponse) -> None:
        self._response = response

    def get(
        self, url: str, timeout: float | tuple[float, float], **kwargs: Any
    ) -> DummyResponse:
        return self._response

    def post(
        self, url: str, timeout: float | tuple[float, float], **kwargs: Any
    ) -> DummyResponse:
        return self._response


def test_do_request_success() -> None:
    """Successful call returns parsed JSON and empty error string."""
    session = DummySession(DummyResponse(200, text="{}", json_data={"a": 1}))
    data, err = pl._do_request(
        cast(requests.Session, session), "http://example.org", delay=0
    )
    assert data == {"a": 1}
    assert err == ""


def test_do_request_404() -> None:
    """404 response is reported as 'PMID not found'."""
    session = DummySession(DummyResponse(404, text="not found", json_data={}))
    data, err = pl._do_request(
        cast(requests.Session, session), "http://example.org", delay=0
    )
    assert data is None
    assert err == "PMID not found"


def test_fetch_pubmed_uses_cfg(monkeypatch) -> None:
    cfg = PubMedCfg(
        base="https://example.org/eutils",
        timeout_connect=1,
        timeout_read=2,
        retries=4,
        encodings=["utf-8"],
    )
    captured: dict[str, Any] = {}

    def fake_do_request(
        session,
        url,
        sleep,
        expect_json=True,
        retries=2,
        method="GET",
        timeout=10,
        **kwargs,
    ):
        captured.update(
            {"url": url, "sleep": sleep, "retries": retries, "timeout": timeout}
        )
        return "<root></root>", ""

    sleeps: list[float] = []
    monkeypatch.setattr(pl, "_do_request", fake_do_request)
    monkeypatch.setattr(pl, "sleep", lambda s: sleeps.append(s))

    session = requests.Session()
    res = pl.fetch_pubmed(session, "1", 0.5, cfg=cfg)
    assert res["PubMed.Error"] == "No PubmedArticle"
    assert (
        captured["url"]
        == "https://example.org/eutils/efetch.fcgi?db=pubmed&id=1&retmode=xml"
    )
    assert captured["timeout"] == (1, 2)
    assert captured["retries"] == 4
    assert captured["sleep"] == 0.5
    assert sleeps == []


def test_fetch_openalex_uses_cfg(monkeypatch) -> None:
    cfg = Config(
        openalex=OpenAlexCfg(
            base="https://example.org",
            timeout_connect=1,
            timeout_read=2,
            retries=4,
            rps=10,
            mailto="info@example.org",
        )
    )

    captured: dict[str, Any] = {}

    def fake_do_request(
        session,
        url,
        sleep,
        expect_json=True,
        retries=2,
        method="GET",
        timeout=10,
        **kwargs,
    ):
        captured.update(
            {"url": url, "sleep": sleep, "retries": retries, "timeout": timeout}
        )
        return {
            "type": "article",
            "type_crossref": "journal-article",
            "genre": "journal",
            "id": "123",
            "host_venue": {"display_name": "Venue"},
            "mesh": [],
        }, ""

    monkeypatch.setattr(pl, "_do_request", fake_do_request)

    class DummyLimiter(rl.RateLimiter):
        def __init__(self) -> None:
            super().__init__(rps=1)
            self.calls = 0

        def acquire(self) -> None:  # type: ignore[override]
            self.calls += 1

    limiter = DummyLimiter()
    session = requests.Session()
    res = pl.fetch_openalex(session, "1", cfg=cfg.openalex, limiter=limiter)
    assert res["OpenAlex.Error"] == ""
    assert (
        captured["url"] == "https://example.org/works/pmid:1?mailto=info%40example.org"
    )
    assert captured["timeout"] == (1, 2)
    assert captured["retries"] == 4
    assert captured["sleep"] == pytest.approx(0.1)
    assert limiter.calls == 1


def test_fetch_crossref_uses_cfg(monkeypatch) -> None:
    cfg = Config(
        crossref=CrossRefCfg(
            base="https://api.example.org",
            timeout_connect=2,
            timeout_read=3,
            retries=5,
            rps=8,
            mailto="info@example.org",
        )
    )

    captured: dict[str, Any] = {}

    def fake_do_request(
        session,
        url,
        sleep,
        expect_json=True,
        retries=2,
        method="GET",
        timeout=10,
        **kwargs,
    ):
        captured.update(
            {"url": url, "sleep": sleep, "retries": retries, "timeout": timeout}
        )
        return {"message": {}}, ""

    monkeypatch.setattr(pl, "_do_request", fake_do_request)

    class DummyLimiter(rl.RateLimiter):
        def __init__(self) -> None:
            super().__init__(rps=1)
            self.calls = 0

        def acquire(self) -> None:  # type: ignore[override]
            self.calls += 1

    limiter = DummyLimiter()
    session = requests.Session()
    res = pl.fetch_crossref(session, "10.1000/xyz", cfg=cfg.crossref, limiter=limiter)
    assert res["crossref.Error"] == ""
    assert (
        captured["url"]
        == "https://api.example.org/works/10.1000%2Fxyz?mailto=info%40example.org"
    )
    assert captured["timeout"] == (2, 3)
    assert captured["retries"] == 5
    assert captured["sleep"] == pytest.approx(1 / 8)
    assert limiter.calls == 1


def test_fetch_semantic_scholar_uses_cfg(monkeypatch) -> None:
    cfg = SemanticScholarCfg(
        base="https://api.example.org/v1",
        timeout_connect=1,
        timeout_read=2,
        retries=4,
        encodings=["utf-8"],
    )
    captured: dict[str, Any] = {}

    def fake_do_request(
        session,
        url,
        sleep,
        expect_json=True,
        retries=2,
        method="GET",
        timeout=10,
        **kwargs,
    ):
        captured.update(
            {"url": url, "sleep": sleep, "retries": retries, "timeout": timeout}
        )
        return {"publicationTypes": [], "externalIds": {}}, ""

    sleeps: list[float] = []
    monkeypatch.setattr(pl, "_do_request", fake_do_request)
    monkeypatch.setattr(pl, "sleep", lambda s: sleeps.append(s))

    session = requests.Session()
    res = pl.fetch_semantic_scholar(session, "1", 0.2, cfg=cfg)
    assert res["scholar.Error"] == ""
    assert captured["url"] == "https://api.example.org/v1/paper/PMID:1"
    assert captured["timeout"] == (1, 2)
    assert captured["retries"] == 4
    assert captured["sleep"] == pytest.approx(1.0)
    assert sleeps == []
