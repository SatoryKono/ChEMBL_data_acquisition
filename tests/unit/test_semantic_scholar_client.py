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
        captured["delay"] = delay
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
    cfg = SemanticScholarCfg(api_key="  secret-token  ", rps=2, delay=3.0)

    limiter = DummyLimiter()

    result = semantic_scholar.fetch_semantic_scholar(
        session,
        "123456",
        cfg=cfg,
        limiter=limiter,
    )

    assert result["scholar.Error"] == ""
    headers = captured.get("headers")
    assert isinstance(headers, dict)
    assert headers["Accept"] == "application/json"
    assert headers["x-api-key"] == "secret-token"
    assert limiter.calls == 1
    delay = captured.get("delay")
    assert delay == pytest.approx(0.5)


@pytest.mark.unit
def test_fetch_semantic_scholar_batch__omits_api_key_header_when_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def _fake_do_request(session, url, delay, *, headers, **kwargs):  # type: ignore[override]
        captured["headers"] = headers
        captured["delay"] = delay
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
    cfg = SemanticScholarCfg(delay=2.5)

    results = semantic_scholar.fetch_semantic_scholar_batch(
        session, ["123456"], cfg=cfg, limiter=limiter
    )

    assert results
    assert results[0]["scholar.Error"] == ""
    headers = captured.get("headers")
    assert isinstance(headers, dict)
    assert headers["Accept"] == "application/json"
    assert "x-api-key" not in headers
    delay = captured.get("delay")
    assert delay == pytest.approx(2.5)


@pytest.mark.unit
@pytest.mark.parametrize(
    "cfg_kwargs,expected",
    [
        ({"rps": 4, "delay": 3.0}, 0.25),
        ({"rps": None, "delay": 1.7}, 1.7),
    ],
)
def test_fetch_semantic_scholar_batch__delay_resolution(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
    cfg_kwargs: dict[str, object],
    expected: float,
) -> None:
    captured: dict[str, float] = {}

    def _fake_do_request(session, url, delay, **kwargs):  # type: ignore[override]
        captured["delay"] = delay
        return ([{"pmid": "1"}], "")

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    session = requests.Session()
    cfg = SemanticScholarCfg(**cfg_kwargs)

    semantic_scholar.fetch_semantic_scholar_batch(session, ["1"], cfg=cfg)

    assert captured["delay"] == pytest.approx(expected)

    def _error_do_request(session, url, delay, **kwargs):  # type: ignore[override]
        return None, "HTTP 500 Internal Server Error"

    monkeypatch.setattr(semantic_scholar, "_do_request", _error_do_request)

    session = requests.Session()

    pmids = ["123456", "789012"]
    with caplog.at_level("WARNING"):
        results = semantic_scholar.fetch_semantic_scholar_batch(session, pmids, 0.0)

    assert all(
        entry["scholar.Error"] == "HTTP 500 Internal Server Error" for entry in results
    )
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


@pytest.mark.unit
def test_fetch_semantic_scholar_batch__access_denied_logs_once(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    def _fake_do_request(session, url, delay, **kwargs):  # type: ignore[override]
        return None, 'HTTP 403: {"message":"Forbidden"}'

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    session = requests.Session()
    pmids = ["1319495", "15317460", "16078848"]

    with caplog.at_level("WARNING"):
        results = semantic_scholar.fetch_semantic_scholar_batch(session, pmids, 0.0)

    assert [entry["scholar.PMID"] for entry in results] == pmids
    for entry in results:
        assert "Semantic Scholar API access denied" in entry["scholar.Error"]
        assert "HTTP 403" in entry["scholar.Error"]

    access_logs = [
        record for record in caplog.records if "access denial" in record.getMessage()
    ]
    assert len(access_logs) == 1
    assert str(len(pmids)) in access_logs[0].getMessage()


@pytest.mark.unit
def test_fetch_semantic_scholar__access_denied_hint(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    def _fake_do_request(session, url, delay, **kwargs):  # type: ignore[override]
        return None, 'HTTP 403: {"message":"Forbidden"}'

    monkeypatch.setattr(semantic_scholar, "_do_request", _fake_do_request)

    session = requests.Session()

    with caplog.at_level("WARNING"):
        result = semantic_scholar.fetch_semantic_scholar(session, "17827018", 0.0)

    assert result["scholar.PMID"] == "17827018"
    assert "Semantic Scholar API access denied" in result["scholar.Error"]
    assert "HTTP 403" in result["scholar.Error"]
    assert any(
        "access denial" in record.getMessage() for record in caplog.records
    )


@pytest.mark.unit
@pytest.mark.parametrize(
    "message,expected",
    [
        ("HTTP 403: Forbidden", True),
        ("http 401: unauthorized", True),
        ("Request failed: Forbidden", True),
        ("HTTP 500: Internal Server Error", False),
        ("", False),
        (None, False),
    ],
)
def test_is_access_denied_error__cases(
    message: str | None,
    expected: bool,
) -> None:
    assert semantic_scholar.is_access_denied_error(message) is expected
