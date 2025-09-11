"""Tests for :mod:`library.chembl_client`."""

from __future__ import annotations

from typing import Any

import random
import requests
import responses  # type: ignore[import-not-found]
import time


import pytest


from cachetools import TTLCache  # type: ignore[import-untyped]


from library import chembl_client



from library.chembl_client import clear_cache, init_session, request_json
from library.config import ApiCfg, ChemblCfg, RetryCfg
import library.rate_limiter as rl


class DummyResponse:
    def __enter__(self) -> "DummyResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        pass

    def raise_for_status(self) -> None:  # pragma: no cover - no error
        pass

    def json(self) -> dict[str, Any]:
        return {"ok": True}


class DummySession:
    def __init__(self, *, failures: int = 0) -> None:
        self.timeout: Any = None
        self.calls: list[int] = []
        self.failures = failures

    def get(self, url: str, timeout: Any) -> DummyResponse:
        self.calls.append(id(self))
        self.timeout = timeout
        if len(self.calls) <= self.failures:
            raise requests.RequestException("boom")
        return DummyResponse()


def test_init_session_sets_user_agent(monkeypatch) -> None:
    """Session should include the configured ``User-Agent`` header."""
    monkeypatch.setattr("library.chembl_client._session", None)
    cfg = ApiCfg(user_agent="test-agent/1.0 (mailto:test@example.org)")
    init_session(cfg, RetryCfg())
    session = chembl_client._session
    assert session is not None
    assert session.headers.get("User-Agent") == cfg.user_agent


class FakeTime:
    def __init__(self) -> None:
        self.now = 0.0
        self.sleeps: list[float] = []

    def monotonic(self) -> float:
        return self.now

    def sleep(self, delay: float) -> None:
        self.sleeps.append(delay)
        self.now += delay


@responses.activate
def test_request_json_uses_cfg(monkeypatch) -> None:
    session = DummySession(failures=1)
    monkeypatch.setattr("library.chembl_client._session", session)
    clear_cache()

    sleep_times: list[float] = []
    monkeypatch.setattr(time, "sleep", lambda t: sleep_times.append(t))
    monkeypatch.setattr(random, "uniform", lambda a, b: 0.0)

    def fail_session(*args, **kwargs):  # pragma: no cover - ensure no new session
        raise AssertionError("requests.Session should not be called")

    monkeypatch.setattr(requests, "Session", fail_session)

    cfg = ApiCfg(timeout_connect=1, timeout_read=2, retries=2, backoff_factor=0.5)
    request_json("http://example.com", cfg=cfg)

    assert session.timeout == (1, 2)
    assert len(session.calls) == 2
    assert sleep_times == [0.5]


def test_request_json_backoff_grows(monkeypatch) -> None:
    session = DummySession(failures=2)
    monkeypatch.setattr("library.chembl_client._session", session)
    clear_cache()

    sleep_times: list[float] = []
    monkeypatch.setattr(time, "sleep", lambda t: sleep_times.append(t))
    monkeypatch.setattr(random, "uniform", lambda a, b: 0.0)

    cfg = ApiCfg(timeout_connect=1, timeout_read=2, retries=3, backoff_factor=1)
    request_json("http://example.com", cfg=cfg)

    assert sleep_times == [1.0, 2.0]


def test_request_json_reuses_session(monkeypatch) -> None:
    dummy = DummySession()

    def fake_session_with_retry(api: ApiCfg, retry: RetryCfg) -> DummySession:
        return dummy

    monkeypatch.setattr(
        "library.chembl_client.session_with_retry", fake_session_with_retry
    )
    init_session(ApiCfg(), RetryCfg())

    def fail_session(*args, **kwargs):  # pragma: no cover - ensure no new session
        raise AssertionError("requests.Session should not be called")

    monkeypatch.setattr(requests, "Session", fail_session)
    clear_cache()

    request_json("http://example.com/1", cfg=ApiCfg())
    request_json("http://example.com/2", cfg=ApiCfg())

    assert dummy.calls == [id(dummy), id(dummy)]


@responses.activate
def test_request_json_cache(monkeypatch) -> None:
    clear_cache()
    monkeypatch.setattr("library.chembl_client._session", None)
    url = "http://example.com/cache"
    responses.add(responses.GET, url, json={"ok": True}, status=200)

    request_json(url, cfg=ApiCfg())
    assert url in chembl_client._CACHE

    request_json(url, cfg=ApiCfg())

    assert url in chembl_client._CACHE

    assert len(responses.calls) == 1


@responses.activate
def test_request_json_cache_ttl_expiration(monkeypatch) -> None:
    timer = [0.0]

    cache: TTLCache[str, Any] = TTLCache(maxsize=2, ttl=1, timer=lambda: timer[0])
    monkeypatch.setattr(chembl_client, "_CACHE", cache)

    monkeypatch.setattr("library.chembl_client._session", None)
    chembl_cfg = ChemblCfg(cache_ttl=1)
    init_session(ApiCfg(), RetryCfg(), chembl_cfg)
    assert chembl_client._CACHE.ttl == chembl_cfg.cache_ttl
    chembl_client._CACHE = TTLCache(
        maxsize=1024, ttl=chembl_cfg.cache_ttl, timer=lambda: timer[0]
    )
    clear_cache()
    url = "http://example.com/ttl"
    responses.add(responses.GET, url, json={"ok": True}, status=200)

    request_json(url, cfg=ApiCfg())

    # Advance time beyond the TTL to force expiration of the cached entry.
    timer[0] = 2.0
    responses.add(responses.GET, url, json={"ok": True}, status=200)
    request_json(url, cfg=ApiCfg())

    # Two HTTP calls should have occurred because the cache entry expired.
    assert len(responses.calls) == 2


@responses.activate
def test_request_json_preserves_original_error_message(monkeypatch) -> None:
    """Ensure the raised error retains status code and URL."""

    clear_cache()
    monkeypatch.setattr("library.chembl_client._session", None)
    url = "http://example.com/notfound"
    responses.add(responses.GET, url, status=404)

    cfg = ApiCfg(retries=1)
    with pytest.raises(requests.HTTPError) as exc_info:
        request_json(url, cfg=cfg)

    message = str(exc_info.value)
    assert "404" in message
    assert url in message


def test_clear_cache(monkeypatch) -> None:
    cache: TTLCache[str, Any] = TTLCache(maxsize=2, ttl=100)
    monkeypatch.setattr(chembl_client, "_CACHE", cache)
    chembl_client._CACHE["x"] = {"ok": True}
    clear_cache()
    assert len(chembl_client._CACHE) == 0


def test_request_json_rate_limiter_blocks(monkeypatch) -> None:
    fake_time = FakeTime()
    monkeypatch.setattr(rl, "time", fake_time)
    monkeypatch.setattr(rl, "sleep", fake_time.sleep)
    with rl._limiters_lock:
        rl._limiters.clear()
    session = DummySession()
    monkeypatch.setattr("library.chembl_client._session", session)
    clear_cache()
    cfg = ApiCfg(rps=1, burst=1)
    request_json("http://example.com/1", cfg=cfg)
    request_json("http://example.com/2", cfg=cfg)
    assert fake_time.sleeps == [1.0]
    with rl._limiters_lock:
        rl._limiters.clear()

