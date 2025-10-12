from __future__ import annotations

import pytest
import requests

from library.clients import semantic_scholar
from library.config.models import SemanticScholarCfg


@pytest.mark.unit
def test_fetch_semantic_scholar__includes_api_key_header(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def _fake_do_request(session, url, delay, *, headers, **kwargs):  # type: ignore[override]
        captured["headers"] = headers
        return (
            {
                "paperId": "S2:123",
                "externalIds": {"DOI": "10.1234/example"},
                "publicationTypes": ["Journal Article"],
                "venue": "Example Journal",
            },
            "",
        )

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    class DummyLimiter:
        def __init__(self) -> None:
            self.calls = 0

        def acquire(self) -> None:
            self.calls += 1

    session = requests.Session()
    cfg = SemanticScholarCfg(api_key="  secret-token  ")

    limiter = DummyLimiter()

    result = semantic_scholar.fetch_semantic_scholar(
        session,
        "123456",
        0.0,
        cfg=cfg,
        limiter=limiter,
    )

    assert result["scholar.Error"] == ""
    headers = captured.get("headers")
    assert isinstance(headers, dict)
    assert headers["Accept"] == "application/json"
    assert headers["x-api-key"] == "secret-token"
    assert limiter.calls == 1


@pytest.mark.unit
def test_fetch_semantic_scholar__logs_warning_on_error(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    def _fake_do_request(session, url, delay, *, headers, **kwargs):  # type: ignore[override]
        return ({}, "HTTP 429 Too Many Requests")

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    session = requests.Session()

    with caplog.at_level("WARNING"):
        result = semantic_scholar.fetch_semantic_scholar(session, "123456", 0.0)

    assert result["scholar.Error"] == "HTTP 429 Too Many Requests"
    assert any(
        "123456" in record.getMessage() and "429" in record.getMessage()
        for record in caplog.records
    )


@pytest.mark.unit
def test_fetch_semantic_scholar__logs_warning_on_invalid_response(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    def _fake_do_request(session, url, delay, *, headers, **kwargs):  # type: ignore[override]
        return ("unexpected", "")

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    session = requests.Session()

    with caplog.at_level("WARNING"):
        result = semantic_scholar.fetch_semantic_scholar(session, "123456", 0.0)

    assert result["scholar.Error"] == "Invalid response"
    assert any(
        "123456" in record.getMessage() and "invalid" in record.getMessage().lower()
        for record in caplog.records
    )


@pytest.mark.unit
def test_fetch_semantic_scholar_batch__omits_api_key_header_when_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def _fake_do_request(session, url, delay, *, headers, **kwargs):  # type: ignore[override]
        captured["headers"] = headers
        return (
            [
                {
                    "paperId": "S2:123",
                    "externalIds": {"DOI": "10.1234/example"},
                    "publicationTypes": ["Journal Article"],
                    "venue": "Example Journal",
                    "pmid": "123456",
                }
            ],
            "",
        )

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    class DummyLimiter:
        def __init__(self) -> None:
            self.calls = 0

        def acquire(self) -> None:
            self.calls += 1

    session = requests.Session()
    limiter = DummyLimiter()

    results = semantic_scholar.fetch_semantic_scholar_batch(
        session,
        ["123456"],
        0.0,
        limiter=limiter,
    )

    assert results
    assert results[0]["scholar.Error"] == ""
    headers = captured.get("headers")
    assert isinstance(headers, dict)
    assert headers["Accept"] == "application/json"
    assert "x-api-key" not in headers
    assert limiter.calls == 1


@pytest.mark.unit
def test_fetch_semantic_scholar__skips_limiter_when_rps_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _fake_do_request(session, url, delay, **kwargs):  # type: ignore[override]
        return (
            {
                "paperId": "S2:123",
                "externalIds": {"DOI": "10.1234/example"},
                "publicationTypes": ["Journal Article"],
                "venue": "Example Journal",
            },
            "",
        )

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    def _fail_get_limiter(*args, **kwargs):  # type: ignore[no-untyped-def]
        raise AssertionError("get_limiter should not be called when rps is None")

    monkeypatch.setattr(semantic_scholar, "get_limiter", _fail_get_limiter)

    session = requests.Session()
    cfg = SemanticScholarCfg()

    result = semantic_scholar.fetch_semantic_scholar(session, "123456", 0.0, cfg=cfg)

    assert result["scholar.Error"] == ""


@pytest.mark.unit
def test_fetch_semantic_scholar__creates_limiter_from_config(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _fake_do_request(session, url, delay, **kwargs):  # type: ignore[override]
        return (
            {
                "paperId": "S2:123",
                "externalIds": {"DOI": "10.1234/example"},
                "publicationTypes": ["Journal Article"],
                "venue": "Example Journal",
            },
            "",
        )

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    class DummyLimiter:
        def __init__(self) -> None:
            self.calls = 0

        def acquire(self) -> None:
            self.calls += 1

    limiter = DummyLimiter()

    def _fake_get_limiter(name, rps, burst):  # type: ignore[no-untyped-def]
        assert name == "semantic_scholar"
        assert rps == 3
        assert burst == 4
        return limiter

    monkeypatch.setattr(semantic_scholar, "get_limiter", _fake_get_limiter)

    session = requests.Session()
    cfg = SemanticScholarCfg(rps=3, burst=4)

    result = semantic_scholar.fetch_semantic_scholar(session, "123456", 0.0, cfg=cfg)

    assert result["scholar.Error"] == ""
    assert limiter.calls == 1


@pytest.mark.unit
def test_fetch_semantic_scholar_batch__logs_warning_on_error(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    def _fake_do_request(session, url, delay, *, headers, **kwargs):  # type: ignore[override]
        return ([], "HTTP 500 Internal Server Error")

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    session = requests.Session()

    pmids = ["123456", "789012"]
    with caplog.at_level("WARNING"):
        results = semantic_scholar.fetch_semantic_scholar_batch(session, pmids, 0.0)

    assert all(entry["scholar.Error"] == "HTTP 500 Internal Server Error" for entry in results)
    for pmid in pmids:
        assert any(
            pmid in record.getMessage() and "500" in record.getMessage()
            for record in caplog.records
        )


@pytest.mark.unit
def test_fetch_semantic_scholar_batch__logs_warning_on_invalid_response(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    def _fake_do_request(session, url, delay, *, headers, **kwargs):  # type: ignore[override]
        return ({"unexpected": "structure"}, "")

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    session = requests.Session()
    pmids = ["123456"]

    with caplog.at_level("WARNING"):
        results = semantic_scholar.fetch_semantic_scholar_batch(session, pmids, 0.0)

    assert results[0]["scholar.Error"] == "Invalid batch response format"
    assert any(
        "123456" in record.getMessage() and "invalid" in record.getMessage().lower()
        for record in caplog.records
    )


@pytest.mark.unit
def test_semantic_scholar_cfg__rejects_blank_api_key() -> None:
    with pytest.raises(ValueError):
        SemanticScholarCfg(api_key="   ")
