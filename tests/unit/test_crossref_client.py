from __future__ import annotations

from typing import Any
from urllib.parse import quote

import pytest
import requests

from library.clients import crossref
from library.common.log import logger
from library.config.models import CrossRefCfg, RetryCfg


class _LimiterStub:
    def __init__(self) -> None:
        self.acquire_calls = 0

    def acquire(self) -> None:
        self.acquire_calls += 1


def _log_request_start(url: str) -> None:
    """Emit a deterministic ``request_start`` log entry for tests."""

    logger.info(
        "request_start",
        extra={"stage": "request_start", "url": url, "attempt": 1},
    )


@pytest.fixture
def limiter_stub(monkeypatch: pytest.MonkeyPatch) -> tuple[_LimiterStub, list[tuple[str, int, int]]]:
    limiter = _LimiterStub()
    calls: list[tuple[str, int, int]] = []

    def _fake_get_limiter(name: str, rps: int, burst: int) -> _LimiterStub:
        calls.append((name, rps, burst))
        return limiter

    monkeypatch.setattr(crossref, "get_limiter", _fake_get_limiter)
    return limiter, calls


@pytest.mark.unit
def test_fetch_crossref__success(
    limiter_stub: tuple[_LimiterStub, list[tuple[str, int, int]]],
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    limiter, limiter_calls = limiter_stub
    session = requests.Session()
    cfg = CrossRefCfg(rps=6, burst=3, retries=2)
    doi = "10.1234/example doi"
    captured: dict[str, Any] = {}

    def _fake_do_request(
        session_arg: requests.Session,
        url: str,
        delay: float,
        *,
        retries: int,
        timeout: tuple[float, float],
        retry_cfg: RetryCfg | None,
    ) -> tuple[dict[str, Any], str]:
        _log_request_start(url)
        captured.update(
            {
                "session": session_arg,
                "url": url,
                "delay": delay,
                "retries": retries,
                "timeout": timeout,
                "retry_cfg": retry_cfg,
            }
        )
        return {"message": "payload"}, ""

    monkeypatch.setattr(crossref, "_do_request", _fake_do_request)

    with caplog.at_level("INFO", logger="chembl"):
        data, error = crossref.fetch_crossref(session, doi, cfg=cfg)

    assert data == {"message": "payload"}
    assert error == ""
    assert limiter.acquire_calls == 1
    assert limiter_calls == [("crossref", cfg.rps, cfg.burst)]
    assert captured["session"] is session
    expected_url = (
        f"{cfg.base.rstrip('/')}/works/{quote(doi, safe='')}?mailto={quote(cfg.mailto)}"
    )
    assert captured["url"] == expected_url
    assert captured["delay"] == pytest.approx(1 / cfg.rps)
    assert captured["retries"] == cfg.retries
    assert captured["timeout"] == (cfg.timeout_connect, cfg.timeout_read)
    assert captured["retry_cfg"] is None
    assert any("request_start" in record.getMessage() for record in caplog.records)
    assert any("request_ok" in record.getMessage() for record in caplog.records)


@pytest.mark.unit
def test_fetch_crossref__returns_404_error(
    limiter_stub: tuple[_LimiterStub, list[tuple[str, int, int]]],
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    limiter, _ = limiter_stub
    session = requests.Session()
    cfg = CrossRefCfg()
    call_count = 0

    def _fake_do_request(*args: Any, **kwargs: Any) -> tuple[None, str]:
        _log_request_start(kwargs.get("url") if "url" in kwargs else args[1])
        nonlocal call_count
        call_count += 1
        return None, "HTTP 404: Not Found"

    monkeypatch.setattr(crossref, "_do_request", _fake_do_request)

    with caplog.at_level("INFO", logger="chembl"):
        data, error = crossref.fetch_crossref(session, "10.9999/not-found", cfg=cfg)

    assert data is None
    assert error == "HTTP 404: Not Found"
    assert limiter.acquire_calls == 1
    assert call_count == 1
    assert any("request_start" in record.getMessage() for record in caplog.records)
    assert any("request_fail" in record.getMessage() for record in caplog.records)


@pytest.mark.unit
def test_fetch_crossref__retry_after_5xx(
    limiter_stub: tuple[_LimiterStub, list[tuple[str, int, int]]],
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    limiter, _ = limiter_stub
    session = requests.Session()
    cfg = CrossRefCfg(rps=8, retries=1)
    retry_cfg = RetryCfg(max_attempts=2, backoff_factor=0.5, backoff_cap=2.0)
    call_args: dict[str, Any] = {}
    call_count = 0

    def _fake_do_request(
        session_arg: requests.Session,
        url: str,
        delay: float,
        *,
        retries: int,
        timeout: tuple[float, float],
        retry_cfg: RetryCfg | None,
    ) -> tuple[None, str]:
        _log_request_start(url)
        nonlocal call_count
        call_count += 1
        call_args.update(
            {
                "session": session_arg,
                "url": url,
                "delay": delay,
                "retries": retries,
                "timeout": timeout,
                "retry_cfg": retry_cfg,
            }
        )
        return None, "HTTP 503: Service Unavailable"

    monkeypatch.setattr(crossref, "_do_request", _fake_do_request)

    with caplog.at_level("INFO", logger="chembl"):
        data, error = crossref.fetch_crossref(
            session, "10.1234/retry", cfg=cfg, retry_cfg=retry_cfg
        )

    assert data is None
    assert error == "HTTP 503: Service Unavailable"
    assert limiter.acquire_calls == 1
    assert call_count == 1
    assert call_args["session"] is session
    assert call_args["delay"] == pytest.approx(1 / cfg.rps)
    assert call_args["retries"] == cfg.retries
    assert call_args["timeout"] == (cfg.timeout_connect, cfg.timeout_read)
    assert call_args["retry_cfg"] is retry_cfg
    assert any("request_start" in record.getMessage() for record in caplog.records)
    assert any("request_fail" in record.getMessage() for record in caplog.records)


@pytest.mark.unit
def test_fetch_crossref__timeout_error(
    limiter_stub: tuple[_LimiterStub, list[tuple[str, int, int]]],
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    limiter, _ = limiter_stub
    session = requests.Session()
    cfg = CrossRefCfg()
    call_count = 0

    def _fake_do_request(*args: Any, **kwargs: Any) -> tuple[None, str]:
        _log_request_start(kwargs.get("url") if "url" in kwargs else args[1])
        nonlocal call_count
        call_count += 1
        return None, "Read timed out"

    monkeypatch.setattr(crossref, "_do_request", _fake_do_request)

    with caplog.at_level("INFO", logger="chembl"):
        data, error = crossref.fetch_crossref(session, "10.3210/timeout", cfg=cfg)

    assert data is None
    assert error == "Read timed out"
    assert limiter.acquire_calls == 1
    assert call_count == 1
    assert any("request_start" in record.getMessage() for record in caplog.records)
    assert any("request_fail" in record.getMessage() for record in caplog.records)
