"""Thread-safety tests for :mod:`library.chembl_client`."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Any

import threading

from library.chembl_client import clear_cache, request_json
from library.config import ApiCfg


class DummyResponse:
    def __enter__(self) -> "DummyResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        pass

    def raise_for_status(self) -> None:
        pass

    def json(self) -> dict[str, Any]:
        return {"ok": True}


class DummySession:
    def __init__(self) -> None:
        self.calls: list[int] = []

    def get(self, url: str, timeout: Any) -> DummyResponse:
        self.calls.append(threading.get_ident())
        return DummyResponse()


def test_single_session_created(monkeypatch) -> None:
    clear_cache()
    monkeypatch.setattr("library.chembl_client._session", None)

    dummy = DummySession()
    create_calls = 0

    def fake_session_with_retry(api: ApiCfg, retry: Any) -> DummySession:
        nonlocal create_calls
        create_calls += 1
        return dummy

    monkeypatch.setattr(
        "library.chembl_client.session_with_retry", fake_session_with_retry
    )

    def worker() -> dict[str, Any]:
        return request_json("http://example.com", cfg=ApiCfg())

    with ThreadPoolExecutor(max_workers=5) as pool:
        futures = [pool.submit(worker) for _ in range(5)]
        results = [f.result() for f in futures]

    assert create_calls == 1
    assert all(r == {"ok": True} for r in results)
