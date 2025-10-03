"""Tests for :mod:`library.pubmed.query`."""

from __future__ import annotations

from pathlib import Path
from types import TracebackType
from typing import Any, cast

import pandas as pd
import pytest
import requests

from library.common import rate_limiter as rl

from library.clients import pubmed as pc
from library.config import Config, PubMedCfg, RetryCfg, SemanticScholarCfg

from library.clients import semantic_scholar as ss_client

from library.pubmed import query as pq

DATA_DIR = Path("tests/data")


def test_read_pmids() -> None:
    path = DATA_DIR / "pmids.csv"
    pmids = pq.read_pmids(path)
    expected = pd.DataFrame({"PMID": ["1", "2"]})
    pd.testing.assert_frame_equal(pmids.reset_index(drop=True), expected)


def test_read_pmids_missing_column(tmp_path: Path) -> None:
    """Ensure ``read_pmids`` validates presence of the PMID column."""
    path = tmp_path / "bad.csv"
    path.write_text("ID\n1\n")
    with pytest.raises(ValueError, match="PMID"):
        pq.read_pmids(path)


class DummyResponse:
    def __init__(
        self,
        status: int,
        text: str = "",
        json_data: dict[str, Any] | None = None,
        headers: dict[str, str] | None = None,
    ) -> None:
        self.status_code = status
        self._text = text
        self._json = json_data
        self._headers = headers or {}

    def __enter__(self) -> DummyResponse:  # pragma: no cover - context manager proto
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:  # pragma: no cover
        return None

    @property
    def text(self) -> str:
        return self._text

    def json(self) -> dict[str, Any]:
        if self._json is None:
            raise ValueError("no json")
        return self._json

    @property
    def headers(self) -> dict[str, str]:
        return self._headers


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
    data, err = pc._do_request(
        cast(requests.Session, session), "http://example.org", delay=0
    )
    assert data == {"a": 1}
    assert err == ""


def test_do_request_404() -> None:
    """404 response is reported as 'PMID not found'."""
    session = DummySession(DummyResponse(404, text="not found", json_data={}))
    data, err = pc._do_request(
        cast(requests.Session, session), "http://example.org", delay=0
    )
    assert data is None
    assert err == "PMID not found"


def test_do_request_attempt_count(monkeypatch: pytest.MonkeyPatch) -> None:
    """Retry loop performs ``retries + 1`` total attempts."""

    statuses = [500, 500, 200]
    calls: list[int] = []

    def fake_make_request(
        *args: Any, **kwargs: Any
    ) -> tuple[int, str, Any, str, dict[str, str]]:
        idx = len(calls)
        calls.append(idx)
        status = statuses[idx]
        if status >= 500:
            return status, "error", None, "", {}
        return status, "{}", {"ok": True}, "", {}

    monkeypatch.setattr(pc, "_make_request", fake_make_request)

    data, err = pc._do_request(
        cast(
            requests.Session, DummySession(DummyResponse(200, text="{}", json_data={}))
        ),
        "http://example.org",
        delay=0,
        retries=2,
    )

    assert data == {"ok": True}
    assert err == ""
    assert len(calls) == 3


def test_do_request_retry_after_sleep(monkeypatch: pytest.MonkeyPatch) -> None:
    """Retry loop should respect Retry-After headers between attempts."""

    responses: list[tuple[int, str, Any, str, dict[str, str]]] = [
        (429, "rate limited", None, "", {"Retry-After": "1.5"}),
        (200, "{}", {"ok": True}, "", {}),
    ]

    def fake_make_request(*args: Any, **kwargs: Any) -> tuple[int, str, Any, str, dict[str, str]]:
        return responses.pop(0)

    sleeps: list[float] = []
    monkeypatch.setattr(pc, "_make_request", fake_make_request)
    monkeypatch.setattr(pc, "sleep", lambda value: sleeps.append(value))

    data, err = pc._do_request(
        requests.Session(),
        "http://example.org",
        delay=3.0,
        retries=1,
        retry_cfg=RetryCfg(backoff_factor=2.0),
    )

    assert data == {"ok": True}
    assert err == ""
    assert sleeps == [pytest.approx(1.5)]


def test_do_request_retry_backoff(monkeypatch: pytest.MonkeyPatch) -> None:
    """Retry loop falls back to RetryCfg exponential backoff."""

    responses: list[tuple[int, str, Any, str, dict[str, str]]] = [
        (500, "error", None, "", {}),
        (200, "{}", {"ok": True}, "", {}),
    ]

    def fake_make_request(*args: Any, **kwargs: Any) -> tuple[int, str, Any, str, dict[str, str]]:
        return responses.pop(0)

    sleeps: list[float] = []
    monkeypatch.setattr(pc, "_make_request", fake_make_request)
    monkeypatch.setattr(pc, "sleep", lambda value: sleeps.append(value))

    retry_cfg = RetryCfg(backoff_factor=0.2)
    data, err = pc._do_request(
        requests.Session(),
        "http://example.org",
        delay=0.0,
        retries=1,
        retry_cfg=retry_cfg,
    )

    assert data == {"ok": True}
    assert err == ""
    assert sleeps == [pytest.approx(0.2)]


def test_handle_response_retryable() -> None:
    """500 response triggers a retry on non-final attempts."""
    data, err, retry, retry_after = pc._handle_response(
        "http://example.org", 500, "error", None, "", True, 0, 1, {"Retry-After": "7"}
    )
    assert data is None
    assert err == ""
    assert retry is True
    assert retry_after == pytest.approx(7.0)


def test_handle_response_retryable_last_attempt() -> None:
    """Retryable error returns a message on the last attempt."""
    data, err, retry, retry_after = pc._handle_response(
        "http://example.org", 500, "error", None, "", True, 1, 1, {}
    )
    assert data is None
    assert retry is False
    assert err.startswith("HTTP 500")
    assert retry_after is None


def test_handle_response_parse_error() -> None:
    """Invalid JSON is reported as a failure."""
    data, err, retry, retry_after = pc._handle_response(
        "http://example.org", 200, "", None, "boom", True, 0, 0, {}
    )
    assert data is None
    assert retry is False
    assert err.startswith("Invalid JSON")
    assert retry_after is None


def test_handle_response_success_text() -> None:
    """Text responses are returned when JSON is not expected."""
    data, err, retry, retry_after = pc._handle_response(
        "http://example.org", 200, "hi", "hi", "", False, 0, 0, {}
    )
    assert data == "hi"
    assert err == ""
    assert retry is False
    assert retry_after is None


def test_handle_response_404() -> None:
    """404 error returns a specific message."""
    data, err, retry, retry_after = pc._handle_response(
        "http://example.org", 404, "not found", None, "", True, 0, 0, {}
    )
    assert data is None
    assert err == "PMID not found"
    assert retry is False
    assert retry_after is None


def test_handle_response_400() -> None:
    """400 error is reported as a bad request."""
    data, err, retry, retry_after = pc._handle_response(
        "http://example.org", 400, "bad", None, "", True, 0, 0, {}
    )
    assert data is None
    assert retry is False
    assert err.startswith("Bad request")
    assert retry_after is None


def test_fetch_pubmed_uses_cfg(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = PubMedCfg(
        base="https://example.org/eutils",
        timeout_connect=1,
        timeout_read=2,
        retries=4,
    )
    captured: dict[str, Any] = {}

    def fake_fetch_pubmed_batch(
        session: requests.Session,
        pmids: list[str],
        delay: float,
        cfg: PubMedCfg | None = None,
        *,
        retry_cfg: Any | None = None,
        client: Any | None = None,
    ) -> tuple[str, str]:
        captured.update(
            {
                "url": f"{cfg.base.rstrip('/')}/efetch.fcgi?db=pubmed&id={','.join(pmids)}&retmode=xml",
                "sleep": delay,
                "retries": cfg.retries,
                "timeout": (cfg.timeout_connect, cfg.timeout_read),
            }
        )
        return "<root></root>", ""

    sleeps: list[float] = []
    monkeypatch.setattr(pc, "fetch_pubmed_batch", fake_fetch_pubmed_batch)
    monkeypatch.setattr(pc, "sleep", lambda s: sleeps.append(s))

    session = requests.Session()
    res = pq.fetch_pubmed(session, "1", 0.5, cfg=cfg)
    assert res["PubMed.Error"] == "No PubmedArticle"
    assert (
        captured["url"]
        == "https://example.org/eutils/efetch.fcgi?db=pubmed&id=1&retmode=xml"
    )
    assert captured["timeout"] == (1, 2)
    assert captured["retries"] == 4
    assert captured["sleep"] == 0.5
    assert sleeps == []


def test_fetch_openalex_uses_cfg(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = Config()
    cfg.api.user_agent = "test@example.com"
    cfg.openalex.base = "https://example.org"
    cfg.openalex.timeout_connect = 1
    cfg.openalex.timeout_read = 2
    cfg.openalex.retries = 4
    cfg.openalex.rps = 10
    cfg.openalex.mailto = "info@example.org"

    captured: dict[str, Any] = {}

    def fake_do_request(
        session: requests.Session,
        url: str,
        sleep: float,
        expect_json: bool = True,
        retries: int = 2,
        method: str = "GET",
        timeout: float = 10,
        **kwargs: Any,
    ) -> tuple[dict[str, Any], str]:
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

    monkeypatch.setattr("library.clients.openalex._do_request", fake_do_request)

    class DummyLimiter(rl.RateLimiter):
        def __init__(self) -> None:
            super().__init__(rps=1)
            self.calls = 0

        def acquire(self) -> None:
            self.calls += 1

    limiter = DummyLimiter()
    session = requests.Session()
    res = pq.fetch_openalex(session, "1", cfg=cfg.openalex, limiter=limiter)
    assert res["OpenAlex.Error"] == ""
    assert (
        captured["url"] == "https://example.org/works/pmid:1?mailto=info%40example.org"
    )
    assert captured["timeout"] == (1, 2)
    assert captured["retries"] == 4
    assert captured["sleep"] == pytest.approx(0.1)
    assert limiter.calls == 1


def test_fetch_crossref_uses_cfg(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = Config()
    cfg.api.user_agent = "test@example.com"
    cfg.crossref.base = "https://api.example.org"
    cfg.crossref.timeout_connect = 2
    cfg.crossref.timeout_read = 3
    cfg.crossref.retries = 5
    cfg.crossref.rps = 8
    cfg.crossref.mailto = "info@example.org"

    captured: dict[str, Any] = {}

    def fake_do_request(
        session: requests.Session,
        url: str,
        sleep: float,
        expect_json: bool = True,
        retries: int = 2,
        method: str = "GET",
        timeout: float = 10,
        **kwargs: Any,
    ) -> tuple[dict[str, Any], str]:
        captured.update(
            {"url": url, "sleep": sleep, "retries": retries, "timeout": timeout}
        )
        return {"message": {}}, ""

    monkeypatch.setattr("library.clients.crossref._do_request", fake_do_request)

    class DummyLimiter(rl.RateLimiter):
        def __init__(self) -> None:
            super().__init__(rps=1)
            self.calls = 0

        def acquire(self) -> None:
            self.calls += 1

    limiter = DummyLimiter()
    session = requests.Session()
    res = pq.fetch_crossref(session, "10.1000/xyz", cfg=cfg.crossref, limiter=limiter)
    assert res["crossref.Error"] == ""
    assert (
        captured["url"]
        == "https://api.example.org/works/10.1000%2Fxyz?mailto=info%40example.org"
    )
    assert captured["timeout"] == (2, 3)
    assert captured["retries"] == 5
    assert captured["sleep"] == pytest.approx(1 / 8)
    assert limiter.calls == 1


def test_fetch_semantic_scholar_uses_cfg(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = SemanticScholarCfg(
        base="https://api.example.org/v1",
        timeout_connect=1,
        timeout_read=2,
        retries=4,
    )
    captured: dict[str, Any] = {}

    def fake_do_request(
        session: requests.Session,
        url: str,
        sleep: float,
        expect_json: bool = True,
        retries: int = 2,
        method: str = "GET",
        timeout: float = 10,
        **kwargs: Any,
    ) -> tuple[dict[str, Any], str]:
        captured.update(
            {"url": url, "sleep": sleep, "retries": retries, "timeout": timeout}
        )
        return {"publicationTypes": [], "externalIds": {}}, ""

    sleeps: list[float] = []
    monkeypatch.setattr("library.clients.semantic_scholar._do_request", fake_do_request)
    monkeypatch.setattr(pc, "sleep", lambda s: sleeps.append(s))

    session = requests.Session()
    res = ss_client.fetch_semantic_scholar(session, "1", 0.2, cfg=cfg)
    assert res["scholar.Error"] == ""
    assert captured["url"] == "https://api.example.org/v1/paper/PMID:1"
    assert captured["timeout"] == (1, 2)
    assert captured["retries"] == 4
    assert captured["sleep"] == pytest.approx(1.0)


def test_fetch_semantic_scholar_batch_partial_response(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Batch lookup keeps valid entries when some PMIDs are missing."""

    payload = [
        {
            "externalIds": {"PubMed": "1", "DOI": "10.1000/example"},
            "paperId": "paper-1",
            "publicationTypes": ["Journal"],
            "venue": "Venue",
        },
        None,
        {
            "externalIds": {"PubMed": ["3"]},
            "paperId": "paper-3",
            "publicationTypes": [],
            "venue": "",
        },
    ]

    monkeypatch.setattr(ss_client, "_do_request", lambda *_, **__: (payload, ""))

    session = requests.Session()
    results = ss_client.fetch_semantic_scholar_batch(session, ["1", "2", "3"], 0.0)

    assert results[0]["scholar.Error"] == ""
    assert results[0]["scholar.PMID"] == "1"
    assert results[0]["scholar.SemanticScholarId"] == "paper-1"
    assert results[0]["scholar.DOI"] == "10.1000/example"

    assert results[1]["scholar.Error"] == "Not found"
    assert results[1]["scholar.PMID"] == "2"

    assert results[2]["scholar.Error"] == ""
    assert results[2]["scholar.PMID"] == "3"
    assert results[2]["scholar.SemanticScholarId"] == "paper-3"
