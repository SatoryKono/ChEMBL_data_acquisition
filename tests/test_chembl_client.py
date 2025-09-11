"""Tests for :mod:`library.chembl_client`."""

from __future__ import annotations

from typing import Any

import random
import time

import pytest
import requests
import responses
from cachetools import TTLCache  # type: ignore[import-untyped]

from library.chembl_client import ChemblClient
from library.config import ApiCfg, ChemblCfg, RetryCfg, session_with_retry
import library.rate_limiter as rl


class DummyResponse:
    """Minimal context manager returning a fixed payload."""

    def __enter__(self) -> "DummyResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # pragma: no cover - no-op
        pass

    def raise_for_status(self) -> None:  # pragma: no cover - no error
        pass

    def json(self) -> dict[str, Any]:
        return {"ok": True}


class DummySession:
    """Session recording calls and allowing configured failures."""

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


def test_session_sets_user_agent() -> None:
    """Client session should include the configured ``User-Agent`` header."""
    cfg = ApiCfg(user_agent="test-agent/1.0 (mailto:test@example.org)")
    client = ChemblClient(cfg, RetryCfg())
    assert client.session.headers.get("User-Agent") == cfg.user_agent


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
def test_request_json_uses_cfg(monkeypatch: pytest.MonkeyPatch) -> None:
    session = DummySession(failures=1)
    client = ChemblClient(session=session)

    sleep_times: list[float] = []
    monkeypatch.setattr(time, "sleep", lambda t: sleep_times.append(t))
    monkeypatch.setattr(random, "uniform", lambda a, b: 0.0)

    cfg = ApiCfg(timeout_connect=1, timeout_read=2, retries=2, backoff_factor=0.5)
    client.request_json("http://example.com", cfg=cfg)

    assert session.timeout == (1, 2)
    assert len(session.calls) == 2
    assert sleep_times == [0.5]


def test_request_json_backoff_grows(monkeypatch: pytest.MonkeyPatch) -> None:
    session = DummySession(failures=2)
    client = ChemblClient(session=session)

    sleep_times: list[float] = []
    monkeypatch.setattr(time, "sleep", lambda t: sleep_times.append(t))
    monkeypatch.setattr(random, "uniform", lambda a, b: 0.0)

    cfg = ApiCfg(timeout_connect=1, timeout_read=2, retries=3, backoff_factor=1)
    client.request_json("http://example.com", cfg=cfg)

    assert sleep_times == [1.0, 2.0]


def test_request_json_reuses_session() -> None:
    dummy = DummySession()

    client = ChemblClient(session=dummy)
    client.request_json("http://example.com/1", cfg=ApiCfg())
    client.request_json("http://example.com/2", cfg=ApiCfg())

    assert dummy.calls == [id(dummy), id(dummy)]


@responses.activate
def test_request_json_cache() -> None:
    client = ChemblClient()
    url = "http://example.com/cache"
    responses.add(responses.GET, url, json={"ok": True}, status=200)

    client.request_json(url, cfg=ApiCfg())
    assert url in client.cache

    client.request_json(url, cfg=ApiCfg())
    assert url in client.cache
    assert len(responses.calls) == 1


@responses.activate
def test_request_json_cache_ttl_expiration() -> None:
    timer = [0.0]
    cache: TTLCache[str, dict[str, Any]] = TTLCache(
        maxsize=2, ttl=1, timer=lambda: timer[0]
    )
    chembl_cfg = ChemblCfg(cache_ttl=1)
    client = ChemblClient(chembl=chembl_cfg, cache=cache)

    url = "http://example.com/ttl"
    responses.add(responses.GET, url, json={"ok": True}, status=200)
    client.request_json(url, cfg=ApiCfg())

    timer[0] = 2.0
    responses.add(responses.GET, url, json={"ok": True}, status=200)
    client.request_json(url, cfg=ApiCfg())

    assert len(responses.calls) == 2


@responses.activate
def test_request_json_preserves_original_error_message() -> None:
    client = ChemblClient()
    url = "http://example.com/notfound"
    responses.add(responses.GET, url, status=404)

    cfg = ApiCfg(retries=1)
    with pytest.raises(requests.HTTPError) as exc_info:
        client.request_json(url, cfg=cfg)

    message = str(exc_info.value)
    assert "404" in message
    assert url in message


def test_clear_cache() -> None:
    cache: TTLCache[str, dict[str, Any]] = TTLCache(maxsize=2, ttl=100)
    client = ChemblClient(cache=cache)
    client.cache["x"] = {"ok": True}
    client.clear_cache()
    assert len(client.cache) == 0


def test_request_json_rate_limiter_blocks(monkeypatch: pytest.MonkeyPatch) -> None:
    fake_time = FakeTime()
    monkeypatch.setattr(rl, "time", fake_time)
    monkeypatch.setattr(rl, "sleep", fake_time.sleep)
    with rl._limiters_lock:
        rl._limiters.clear()
    session = DummySession()
    client = ChemblClient(session=session)
    cfg = ApiCfg(rps=1, burst=1)
    client.request_json("http://example.com/1", cfg=cfg)
    client.request_json("http://example.com/2", cfg=cfg)
    assert fake_time.sleeps == [1.0]
    with rl._limiters_lock:
        rl._limiters.clear()


@responses.activate
def test_clients_isolated_caches() -> None:
    url = "http://example.com/data"
    responses.add(responses.GET, url, json={"ok": True}, status=200)
    client1 = ChemblClient()
    client1.request_json(url, cfg=ApiCfg())
    assert len(responses.calls) == 1
    client2 = ChemblClient()
    responses.add(responses.GET, url, json={"ok": True}, status=200)
    client2.request_json(url, cfg=ApiCfg())
    assert len(responses.calls) == 2


@responses.activate
def test_session_with_retry_retries_post() -> None:
    """Ensure POST requests are retried when configured."""

    url = "http://example.com/post"
    responses.add(responses.POST, url, status=500)
    responses.add(responses.POST, url, json={"ok": True}, status=200)

    session = session_with_retry(ApiCfg(), RetryCfg(max_attempts=2, backoff_factor=0))
    response = session.post(url)

    assert response.status_code == 200
    assert len(responses.calls) == 2
