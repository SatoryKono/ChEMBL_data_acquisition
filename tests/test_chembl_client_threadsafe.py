"""Thread-safety tests for :mod:`library.clients.chembl`."""

from __future__ import annotations

import threading
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import pytest

import library.common.rate_limiter as rl
from library.clients import ChemblClient
from library.config import ApiCfg, RetryCfg
from library.clients import chembl as chembl_module

USER_AGENT = "test-agent/1.0 (mailto:test@example.org)"


def api_cfg(**kwargs: Any) -> ApiCfg:
    return ApiCfg(user_agent=USER_AGENT, **kwargs)


class DummyResponse:
    """Minimal response object returning a constant payload."""

    def __enter__(self) -> DummyResponse:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # pragma: no cover - no-op
        pass

    def raise_for_status(self) -> None:  # pragma: no cover - no error
        pass

    def json(self) -> dict[str, Any]:
        return {"ok": True}


class DummySession:
    """HTTP session counting requests and yielding dummy responses."""

    def __init__(self) -> None:
        self.calls = 0

    def get(self, url: str, timeout: Any) -> DummyResponse:
        self.calls += 1
        return DummyResponse()

    def close(self) -> None:  # pragma: no cover - tests supply lightweight sessions
        pass


class ThreadTrackingSession:
    """Session recording the thread identifier for each request."""

    def __init__(self) -> None:
        self.calls: list[int] = []

    def get(self, url: str, timeout: Any) -> DummyResponse:
        self.calls.append(threading.get_ident())
        return DummyResponse()

    def close(self) -> None:  # pragma: no cover - tests supply lightweight sessions
        pass


def test_request_json_threadsafe(monkeypatch) -> None:
    """Concurrent calls should yield the same result as sequential ones."""

    created_sessions: list[DummySession] = []

    def fake_session_with_retry(api: ApiCfg, retry: Any) -> DummySession:
        session = DummySession()
        created_sessions.append(session)
        return session

    monkeypatch.setattr(
        "library.clients.chembl.session_with_retry", fake_session_with_retry
    )
    client = ChemblClient(api_cfg(), RetryCfg())
    client.clear_cache()

    url = "http://example.com/threadsafe"
    cfg = api_cfg()

    # Sequential calls: only the first should trigger a real request.
    sequential = [client.request_json(url, cfg=cfg) for _ in range(5)]
    assert client.session.calls == 1

    client.clear_cache()

    # Parallel calls starting at the same time.
    results: list[dict[str, Any]] = []
    start = threading.Event()

    def worker() -> None:
        start.wait()
        results.append(client.request_json(url, cfg=cfg))

    threads = [threading.Thread(target=worker) for _ in range(5)]
    for thread in threads:
        thread.start()
    start.set()
    for thread in threads:
        thread.join()

    assert results == sequential
    assert len(created_sessions) >= 1


def test_thread_local_sessions_created(monkeypatch) -> None:
    """Each thread should receive its own HTTP session."""

    created_sessions: list[ThreadTrackingSession] = []

    def fake_session_with_retry(api: ApiCfg, retry: Any) -> ThreadTrackingSession:
        session = ThreadTrackingSession()
        created_sessions.append(session)
        return session

    monkeypatch.setattr(
        "library.clients.chembl.session_with_retry", fake_session_with_retry
    )
    client = ChemblClient(api_cfg(), RetryCfg())

    urls = [f"http://example.com/thread/{idx}" for idx in range(5)]

    barrier = threading.Barrier(len(urls))

    def worker(target_url: str) -> dict[str, Any]:
        barrier.wait()
        return client.request_json(target_url, cfg=api_cfg())

    with ThreadPoolExecutor(max_workers=len(urls)) as pool:
        futures = [pool.submit(worker, url) for url in urls]
        results = [future.result() for future in futures]

    assert all(result == {"ok": True} for result in results)
    assert len(created_sessions) == len(urls)
    assert all(len(session.calls) == 1 for session in created_sessions)
    thread_ids = {session.calls[0] for session in created_sessions}
    assert len(thread_ids) == len(urls)


def test_cache_shared_across_threads(monkeypatch) -> None:
    """Cached responses should be reused regardless of the calling thread."""

    created_sessions: list[DummySession] = []

    def fake_session_with_retry(api: ApiCfg, retry: Any) -> DummySession:
        session = DummySession()
        created_sessions.append(session)
        return session

    monkeypatch.setattr(
        "library.clients.chembl.session_with_retry", fake_session_with_retry
    )
    client = ChemblClient(api_cfg(), RetryCfg())
    client.clear_cache()

    url = "http://example.com/data"

    assert client.request_json(url, cfg=api_cfg()) == {"ok": True}
    assert client.session.calls == 1
    assert len(created_sessions) == 1

    def worker() -> dict[str, Any]:
        return client.request_json(url, cfg=api_cfg())

    with ThreadPoolExecutor(max_workers=5) as pool:
        futures = [pool.submit(worker) for _ in range(5)]
        results = [future.result() for future in futures]

    assert all(result == {"ok": True} for result in results)
    assert len(created_sessions) == 1
    assert client.session.calls == 1


def test_concurrent_fetches_no_longer_block(monkeypatch) -> None:
    """Concurrent fetches should proceed in parallel without serialisation."""

    active_lock = threading.Lock()
    active = 0
    max_active = 0

    class SlowResponse:
        def __enter__(self) -> SlowResponse:
            return self

        def __exit__(self, exc_type, exc, tb) -> None:  # pragma: no cover - no-op
            pass

        def raise_for_status(self) -> None:  # pragma: no cover - no error
            pass

        def json(self) -> dict[str, Any]:
            return {"ok": True}

    class SlowSession:
        def get(self, url: str, timeout: Any) -> SlowResponse:
            nonlocal active, max_active
            with active_lock:
                active += 1
                max_active = max(max_active, active)
            try:
                time.sleep(0.05)
                return SlowResponse()
            finally:
                with active_lock:
                    active -= 1

        def close(self) -> None:  # pragma: no cover - tests supply lightweight sessions
            pass

    def fake_session_with_retry(api: ApiCfg, retry: Any) -> SlowSession:
        return SlowSession()

    monkeypatch.setattr(
        "library.clients.chembl.session_with_retry", fake_session_with_retry
    )
    client = ChemblClient(api_cfg(), RetryCfg())
    client.clear_cache()

    urls = [f"http://example.com/stress/{idx}" for idx in range(5)]
    cfg = api_cfg(rps=100, burst=100)
    start = threading.Event()
    results: list[dict[str, Any]] = []

    def worker(target_url: str) -> None:
        start.wait()
        results.append(client.request_json(target_url, cfg=cfg))

    threads = [threading.Thread(target=worker, args=(url,)) for url in urls]
    for thread in threads:
        thread.start()

    start.set()

    for thread in threads:
        thread.join()

    assert len(results) == len(urls)
    assert all(result == {"ok": True} for result in results)
    assert max_active > 1


def test_parallel_clients_respect_global_rate(monkeypatch) -> None:
    """Multiple clients should share the system-wide rate limiter."""

    fake_time = 0.0
    sleep_calls: list[float] = []

    def fake_monotonic() -> float:
        return fake_time

    def fake_sleep(delay: float) -> None:
        nonlocal fake_time
        sleep_calls.append(delay)
        fake_time += delay

    monkeypatch.setattr(rl.time, "monotonic", fake_monotonic)
    monkeypatch.setattr(rl, "sleep", fake_sleep)

    class ImmediateLimiter:
        def acquire(self) -> None:  # pragma: no cover - simple stub
            return None

    monkeypatch.setattr(
        chembl_module,
        "get_limiter",
        lambda *args, **kwargs: ImmediateLimiter(),
    )

    global_limiter = rl.RateLimiter(rps=2, burst=2)

    class RecordingSession:
        def __init__(self) -> None:
            self.calls: list[str] = []

        def get(self, url: str, timeout: Any) -> DummyResponse:
            self.calls.append(url)
            return DummyResponse()

        def close(self) -> None:  # pragma: no cover - compatibility stub
            return None

    cfg = api_cfg(rps=100, burst=100)
    session_one = RecordingSession()
    session_two = RecordingSession()

    client_one = ChemblClient(
        cfg,
        RetryCfg(),
        global_limiter=global_limiter,
        session=session_one,
    )
    client_two = ChemblClient(
        cfg,
        RetryCfg(),
        global_limiter=global_limiter,
        session=session_two,
    )

    urls = [f"http://example.com/global/{idx}" for idx in range(4)]
    barrier = threading.Barrier(len(urls))

    def worker(client: ChemblClient, url: str) -> None:
        barrier.wait()
        client.request_json(url, cfg=cfg)

    threads: list[threading.Thread] = []
    for index, url in enumerate(urls):
        target_client = client_one if index % 2 == 0 else client_two
        thread = threading.Thread(target=worker, args=(target_client, url))
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    assert len(session_one.calls) == 2
    assert len(session_two.calls) == 2
    assert set(session_one.calls) == set(urls[::2])
    assert set(session_two.calls) == set(urls[1::2])
    assert len(sleep_calls) == 2
    assert all(pytest.approx(0.5, rel=1e-6) == call for call in sleep_calls)
    assert fake_time == pytest.approx(1.0, rel=1e-6)
