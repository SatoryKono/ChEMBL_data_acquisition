"""Tests for :mod:`library.chembl_client`."""

from __future__ import annotations

from typing import Any

import random
import requests
import responses
import time
from library.chembl_client import init_session, request_json
from library.config import ApiCfg, RetryCfg


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


@responses.activate
def test_request_json_uses_cfg(monkeypatch) -> None:
    session = DummySession(failures=1)
    monkeypatch.setattr("library.chembl_client._session", session)
    monkeypatch.setattr("library.chembl_client._CACHE", {})

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
    monkeypatch.setattr("library.chembl_client._CACHE", {})

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
    monkeypatch.setattr("library.chembl_client._CACHE", {})

    request_json("http://example.com/1", cfg=ApiCfg())
    request_json("http://example.com/2", cfg=ApiCfg())

    assert dummy.calls == [id(dummy), id(dummy)]
