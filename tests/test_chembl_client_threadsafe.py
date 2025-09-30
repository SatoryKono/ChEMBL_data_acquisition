"""Thread-safety tests for :mod:`library.clients.chembl`."""

from __future__ import annotations

import threading
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any

from library.clients import ChemblClient
from library.config import ApiCfg, RetryCfg

USER_AGENT = "test-agent/1.0 (mailto:test@example.org)"


def api_cfg(**kwargs: Any) -> ApiCfg:
    return ApiCfg(user_agent=USER_AGENT, **kwargs)


class DummyResponse:
    """Minimal response object returning a unique call identifier."""

    def __init__(self, call_no: int) -> None:
        self._call_no = call_no

    def __enter__(self) -> DummyResponse:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # pragma: no cover - no-op
        pass

    def raise_for_status(self) -> None:  # pragma: no cover - no error
        pass

    def json(self) -> dict[str, Any]:
        return {"ok": True}


class DummySession:
    """HTTP session counting requests and yielding incrementing responses."""

    def __init__(self) -> None:
        self.calls = 0

    def get(self, url: str, timeout: Any) -> DummyResponse:
        self.calls += 1
        return DummyResponse(self.calls)


class ThreadTrackingResponse:
    """Response returning a constant payload."""

    def __enter__(self) -> ThreadTrackingResponse:
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


class StressResponse:
    """Response keeping the owning session marked as in use."""

    def __init__(self, session: StressSession, thread_id: int) -> None:
        self._session = session
        self._thread_id = thread_id

    def __enter__(self) -> StressResponse:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # pragma: no cover - no-op
        self._session.release()

    def raise_for_status(self) -> None:  # pragma: no cover - no error
        pass

    def json(self) -> dict[str, Any]:
        time.sleep(0.005)
        return {"thread": self._thread_id}


class StressSession:
    """Session detecting concurrent ``get`` calls from multiple threads."""

    def __init__(self) -> None:
        self.calls: list[int] = []
        self._lock = threading.Lock()
        self._in_use = False

    def release(self) -> None:
        with self._lock:
            self._in_use = False

    def get(self, url: str, timeout: Any) -> StressResponse:
        thread_id = threading.get_ident()
        with self._lock:
            if self._in_use:
                raise RuntimeError("Concurrent access to shared session")
            self._in_use = True
        self.calls.append(thread_id)
        return StressResponse(self, thread_id)


def test_request_json_threadsafe(monkeypatch) -> None:
    """Concurrent calls should yield the same result as sequential ones."""

    session = DummySession()
    client = ChemblClient(api_cfg(), RetryCfg(), session=session)
    client.clear_cache()

    url = "http://example.com/threadsafe"
    cfg = api_cfg()

    # Sequential calls: only the first should trigger a real request.
    sequential = [client.request_json(url, cfg=cfg) for _ in range(5)]
    assert session.calls == 1

    client.clear_cache()
    session.calls = 0

    # Parallel calls starting at the same time.
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


def test_single_session_created(monkeypatch) -> None:
    """Ensure only one HTTP session is initialised across threads."""

    create_calls = 0
    dummy = ThreadTrackingSession()

    def fake_session_with_retry(api: ApiCfg, retry: Any) -> ThreadTrackingSession:
        nonlocal create_calls
        create_calls += 1
        return dummy

    monkeypatch.setattr(
        "library.clients.chembl.session_with_retry", fake_session_with_retry
    )
    client = ChemblClient(api_cfg(), RetryCfg())

    def worker() -> dict[str, Any]:
        return client.request_json("http://example.com", cfg=api_cfg())

    with ThreadPoolExecutor(max_workers=5) as pool:
        futures = [pool.submit(worker) for _ in range(5)]
        results = [f.result() for f in futures]

    assert create_calls == 1
    assert all(r == {"ok": True} for r in results)


def test_cache_shared_across_threads(monkeypatch) -> None:
    client = ChemblClient(api_cfg(), RetryCfg(), session=DummySession())
    client.clear_cache()

    url = "http://example.com/data"
    assert client.request_json(url, cfg=api_cfg()) == {"ok": True}
    assert client.session.calls == 1

    def worker() -> dict[str, Any]:
        return client.request_json(url, cfg=api_cfg())

    with ThreadPoolExecutor(max_workers=5) as pool:
        futures = [pool.submit(worker) for _ in range(5)]
        results = [f.result() for f in futures]

    assert all(r == {"ok": True} for r in results)
    assert client.session.calls == 1


def test_session_lock_prevents_concurrent_access() -> None:
    """High concurrency should serialise access to the shared session."""

    session = StressSession()
    client = ChemblClient(api_cfg(), RetryCfg(), session=session)
    client.clear_cache()

    start = threading.Event()
    urls = [f"http://example.com/stress/{idx}" for idx in range(20)]
    results: list[dict[str, Any]] = []
    errors: list[Exception] = []

    def worker(target_url: str) -> None:
        start.wait()
        try:
            results.append(client.request_json(target_url, cfg=api_cfg()))
        except Exception as exc:  # pragma: no cover - exercised on failure
            errors.append(exc)

    threads = [threading.Thread(target=worker, args=(url,)) for url in urls]
    for thread in threads:
        thread.start()

    start.set()
    for thread in threads:
        thread.join()

    assert not errors
    assert len(results) == len(urls)
    assert session.calls == [result["thread"] for result in results]
    assert len(set(session.calls)) > 1
