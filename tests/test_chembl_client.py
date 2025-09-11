"""Tests for :mod:`library.chembl_client`."""

from __future__ import annotations

from typing import Any

import requests
from library.chembl_client import request_json
from library.config import ApiCfg


class DummyResponse:
    def __enter__(self) -> "DummyResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        pass

    def raise_for_status(self) -> None:  # pragma: no cover - no error
        pass

    def json(self) -> dict[str, Any]:
        return {"ok": True}


class DummyAdapter:
    def __init__(self, max_retries):
        self.retries = max_retries


class DummySession:
    def __init__(self) -> None:
        self.timeout: Any = None
        self.mounted: list[tuple[str, DummyAdapter]] = []

    def __enter__(self) -> "DummySession":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        pass

    def mount(self, prefix: str, adapter: DummyAdapter) -> None:
        self.mounted.append((prefix, adapter))

    def get(self, url: str, timeout: Any) -> DummyResponse:
        self.timeout = timeout
        return DummyResponse()


def test_request_json_uses_cfg(monkeypatch) -> None:
    captured: dict[str, Any] = {}

    def fake_session() -> DummySession:
        sess = DummySession()
        captured["session"] = sess
        return sess

    def fake_adapter(max_retries):
        captured["retries"] = max_retries.total
        captured["backoff"] = max_retries.backoff_factor
        return DummyAdapter(max_retries)

    monkeypatch.setattr(requests, "Session", fake_session)
    monkeypatch.setattr("library.chembl_client.HTTPAdapter", fake_adapter)

    cfg = ApiCfg(timeout_connect=1, timeout_read=2, retries=4, backoff_factor=0.1)
    request_json("http://example.com", cfg=cfg)

    session = captured["session"]
    assert session.timeout == (1, 2)
    assert captured["retries"] == 4
    assert captured["backoff"] == 0.1
