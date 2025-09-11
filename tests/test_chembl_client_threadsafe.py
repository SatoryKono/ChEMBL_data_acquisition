import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import pytest

from library.chembl_client import ChemblClient
from library.config import ApiCfg


class DummyResponse:
    """Minimal response object returning a unique call identifier."""

    def __init__(self, call_no: int) -> None:
        self._call_no = call_no

    def __enter__(self) -> "DummyResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # pragma: no cover - no-op
        pass

    def raise_for_status(self) -> None:  # pragma: no cover - no error
        pass

    def json(self) -> dict[str, Any]:
        return {"call": self._call_no}


class DummySession:
    """HTTP session counting requests and yielding incrementing responses."""

    def __init__(self) -> None:
        self.calls = 0

    def get(self, url: str, timeout: Any) -> DummyResponse:
        self.calls += 1
        return DummyResponse(self.calls)


class ThreadTrackingResponse:
    """Response returning a constant payload."""

    def __enter__(self) -> "ThreadTrackingResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # pragma: no cover - no-op
        pass

    def raise_for_status(self) -> None:  # pragma: no cover - no error
        pass

    def json(self) -> dict[str, Any]:
        return {"ok": True}


class ThreadTrackingSession:
    """Session recording the thread identifier for each request."""

    def __init__(self) -> None:
        self.calls: list[int] = []

    def get(self, url: str, timeout: Any) -> ThreadTrackingResponse:
        self.calls.append(threading.get_ident())
        return ThreadTrackingResponse()


def test_request_json_threadsafe() -> None:
    """Concurrent calls should yield the same result as sequential ones."""

    session = DummySession()
    client = ChemblClient(session=session)

    url = "http://example.com/threadsafe"
    cfg = ApiCfg()

    sequential = [client.request_json(url, cfg=cfg) for _ in range(5)]
    assert session.calls == 1

    client.clear_cache()
    session.calls = 0

    results: list[dict[str, Any]] = []
    start = threading.Event()

    def worker() -> None:
        start.wait()
        results.append(client.request_json(url, cfg=cfg))

    threads = [threading.Thread(target=worker) for _ in range(5)]
    for t in threads:
        t.start()
    start.set()
    for t in threads:
        t.join()

    assert results == sequential


def test_single_session_created(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure only one HTTP session is initialised when constructing a client."""

    dummy = ThreadTrackingSession()
    create_calls = 0

    def fake_session_with_retry(api: ApiCfg, retry: Any) -> ThreadTrackingSession:
        nonlocal create_calls
        create_calls += 1
        return dummy

    monkeypatch.setattr(
        "library.chembl_client.session_with_retry", fake_session_with_retry
    )

    client = ChemblClient()

    def worker() -> dict[str, Any]:
        return client.request_json("http://example.com", cfg=ApiCfg())

    with ThreadPoolExecutor(max_workers=5) as pool:
        futures = [pool.submit(worker) for _ in range(5)]
        results = [f.result() for f in futures]

    assert create_calls == 1
    assert all(r == {"ok": True} for r in results)


def test_cache_shared_across_threads() -> None:
    client = ChemblClient(session=DummySession())
    url = "http://example.com/data"
    assert client.request_json(url, cfg=ApiCfg()) == {"call": 1}

    def worker() -> dict[str, Any]:
        return client.request_json(url, cfg=ApiCfg())

    with ThreadPoolExecutor(max_workers=5) as pool:
        futures = [pool.submit(worker) for _ in range(5)]
        results = [f.result() for f in futures]

    assert all(r == {"call": 1} for r in results)
    assert client.session.calls == 1
