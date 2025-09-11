"""Tests for :mod:`library.chembl_client`."""

from __future__ import annotations

from typing import Any

import requests
import responses
from library.chembl_client import request_json
from library.config import ApiCfg


class CaptureSession(requests.Session):
    def __init__(self) -> None:
        super().__init__()
        self.timeout: Any = None

    def get(self, url: str, timeout: Any | None = None, **kwargs):  # type: ignore[override]
        self.timeout = timeout
        return super().get(url, timeout=timeout, **kwargs)


@responses.activate
def test_request_json_uses_cfg(monkeypatch) -> None:
    captured: dict[str, Any] = {}

    def fake_session() -> CaptureSession:
        sess = CaptureSession()
        captured["session"] = sess
        return sess

    def fake_adapter(max_retries):
        captured["retries"] = max_retries.total
        captured["backoff"] = max_retries.backoff_factor
        return requests.adapters.HTTPAdapter(max_retries=max_retries)

    monkeypatch.setattr(requests, "Session", fake_session)
    monkeypatch.setattr("library.chembl_client.HTTPAdapter", fake_adapter)
    responses.add(responses.GET, "http://example.com", json={"ok": True})

    cfg = ApiCfg(timeout_connect=1, timeout_read=2, retries=4, backoff_factor=0.1)
    data = request_json("http://example.com", cfg=cfg)
    assert data == {"ok": True}

    session = captured["session"]
    assert session.timeout == (1, 2)
    assert captured["retries"] == 4
    assert captured["backoff"] == 0.1
