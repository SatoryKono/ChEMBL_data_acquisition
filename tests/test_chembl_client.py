from __future__ import annotations

import random
import time
from typing import Any
from unittest.mock import MagicMock

import pytest

from library import rate_limiter as rl
from library.chembl_client import ChemblClient
from library.config import ApiCfg, RetryCfg

cachetools = pytest.importorskip("cachetools")
from cachetools import TTLCache  # type: ignore[import-untyped]  # noqa: E402

requests = pytest.importorskip("requests")
responses = pytest.importorskip("responses")

USER_AGENT = "test-agent/1.0 (mailto:test@example.org)"


def api_cfg(**kwargs: Any) -> ApiCfg:
    """Return :class:`ApiCfg` with a default test user agent."""

    return ApiCfg(user_agent=USER_AGENT, **kwargs)


class DummyResponse:
    def __enter__(self) -> DummyResponse:  # pragma: no cover - trivial
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # pragma: no cover - no-op
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


def test_client_sets_user_agent() -> None:
    """Session should include the configured ``User-Agent`` header."""

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
def test_request_json_uses_cfg(monkeypatch) -> None:
    session = DummySession(failures=1)
    client = ChemblClient(api_cfg(), RetryCfg(), session=session)
    client.clear_cache()

    sleep_times: list[float] = []
    monkeypatch.setattr(time, "sleep", lambda t: sleep_times.append(t))
    monkeypatch.setattr(random, "uniform", lambda a, b: 0.0)

    cfg = api_cfg(timeout_connect=1, timeout_read=2, retries=2, backoff_factor=0.5)
    client.request_json("http://example.com", cfg=cfg)

    assert session.timeout == (1, 2)
    assert len(session.calls) == 2
    assert sleep_times == [0.5]


def test_request_json_backoff_grows(monkeypatch) -> None:
    session = DummySession(failures=2)
    client = ChemblClient(api_cfg(), RetryCfg(), session=session)
    client.clear_cache()

    sleep_times: list[float] = []
    monkeypatch.setattr(time, "sleep", lambda t: sleep_times.append(t))
    monkeypatch.setattr(random, "uniform", lambda a, b: 0.0)

    cfg = api_cfg(timeout_connect=1, timeout_read=2, retries=3, backoff_factor=1)
    client.request_json("http://example.com", cfg=cfg)

    assert sleep_times == [1.0, 2.0]


def test_request_json_retries_with_mocked_session() -> None:
    """Failing requests should be retried using the provided session."""

    response = DummyResponse()
    session = MagicMock(spec=requests.Session)
    session.get.side_effect = [requests.RequestException("boom"), response]
    client = ChemblClient(api_cfg(), RetryCfg(), session=session)
    client.clear_cache()

    cfg = api_cfg(retries=2, backoff_factor=0)
    result = client.request_json("http://example.com", cfg=cfg)

    assert result == {"ok": True}
    assert session.get.call_count == 2


def test_request_json_reuses_session(monkeypatch) -> None:
    dummy = DummySession()

    def fake_session_with_retry(api: ApiCfg, retry: RetryCfg) -> DummySession:
        return dummy

    monkeypatch.setattr(
        "library.chembl_client.session_with_retry", fake_session_with_retry
    )
    client = ChemblClient(api_cfg(), RetryCfg())
    client.clear_cache()

    client.request_json("http://example.com/1", cfg=api_cfg())
    client.request_json("http://example.com/2", cfg=api_cfg())

    assert dummy.calls == [id(dummy), id(dummy)]


@responses.activate
def test_request_json_cache(monkeypatch) -> None:
    client = ChemblClient(api_cfg(), RetryCfg())
    client.clear_cache()
    url = "http://example.com/cache"
    responses.add(responses.GET, url, json={"ok": True}, status=200)

    client.request_json(url, cfg=api_cfg())
    assert url in client.cache

    client.request_json(url, cfg=api_cfg())

    assert url in client.cache
    assert len(responses.calls) == 1


@responses.activate
def test_request_json_cache_ttl_expiration(monkeypatch) -> None:
    timer = [0.0]

    cache: TTLCache[str, dict[str, Any]] = TTLCache(
        maxsize=2, ttl=1, timer=lambda: timer[0]
    )

    client = ChemblClient(api_cfg(), RetryCfg())
    client.cache = cache
    client.clear_cache()
    assert client.cache.ttl == 1

    url = "http://example.com/ttl"
    responses.add(responses.GET, url, json={"ok": True}, status=200)

    client.request_json(url, cfg=api_cfg())

    timer[0] = 2.0
    responses.add(responses.GET, url, json={"ok": True}, status=200)
    client.request_json(url, cfg=api_cfg())

    assert len(responses.calls) == 2


@responses.activate
def test_request_json_preserves_original_error_message(monkeypatch) -> None:
    client = ChemblClient(api_cfg(), RetryCfg())
    client.clear_cache()
    url = "http://example.com/notfound"
    responses.add(responses.GET, url, status=404)

    cfg = api_cfg(retries=1)
    with pytest.raises(requests.HTTPError) as exc_info:
        client.request_json(url, cfg=cfg)

    message = str(exc_info.value)
    assert "404" in message
    assert url in message


def test_clear_cache() -> None:
    cache: TTLCache[str, dict[str, Any]] = TTLCache(maxsize=2, ttl=100)
    client = ChemblClient(api_cfg(), RetryCfg())
    client.cache = cache
    client.cache["x"] = {"ok": True}
    client.clear_cache()
    assert len(client.cache) == 0


def test_request_json_rate_limiter_blocks(monkeypatch) -> None:
    fake_time = FakeTime()
    monkeypatch.setattr(rl, "time", fake_time)
    monkeypatch.setattr(rl, "sleep", fake_time.sleep)
    with rl._limiters_lock:
        rl._limiters.clear()
    session = DummySession()
    client = ChemblClient(api_cfg(), RetryCfg(), session=session)
    client.clear_cache()
    cfg = api_cfg(rps=1, burst=1)
    client.request_json("http://example.com/1", cfg=cfg)
    client.request_json("http://example.com/2", cfg=cfg)
    assert fake_time.sleeps == [1.0]
    with rl._limiters_lock:
        rl._limiters.clear()


@responses.activate
def test_clients_do_not_share_cache() -> None:
    url = "http://example.com/data"
    responses.add(responses.GET, url, json={"ok": 1}, status=200)
    client1 = ChemblClient(api_cfg(), RetryCfg())
    client1.request_json(url, cfg=api_cfg())
    responses.add(responses.GET, url, json={"ok": 2}, status=200)
    client2 = ChemblClient(api_cfg(), RetryCfg())
    client2.request_json(url, cfg=api_cfg())
    assert len(responses.calls) == 2
