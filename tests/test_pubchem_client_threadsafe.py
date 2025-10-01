"""Thread-safety tests for :mod:`library.clients.pubchem`."""

from __future__ import annotations

import threading
from concurrent.futures import ThreadPoolExecutor

import pytest

from cachetools import TTLCache

from library import pubchem_library as pl
from library.clients import pubchem as pc
from library.config import ApiCfg, RetryCfg, session_with_retry


class DummyResponse:
    """Simple response object returning a constant payload."""

    status_code = 200

    def __init__(
        self, payload: dict[str, object], close_log: list[object] | None = None
    ) -> None:
        self._payload = payload
        self._close_log = close_log
        self.closed = False

    def __enter__(self) -> DummyResponse:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()
        return None

    def close(self) -> None:
        self.closed = True
        if self._close_log is not None:
            self._close_log.append(True)

    def raise_for_status(self) -> None:  # pragma: no cover - no error
        return None

    def json(self) -> dict[str, object]:
        return self._payload


class NoopLimiter:
    """Limiter stub avoiding waits in tests."""

    def acquire(self) -> None:  # pragma: no cover - simple stub
        return None


class DummySession:
    """Minimal session stub recording close operations."""

    def __init__(self, user_agent: str, closed_log: list[str] | None = None) -> None:
        self.headers = {"User-Agent": user_agent}
        self._closed_log = closed_log
        self.closed = False

    def close(self) -> None:  # pragma: no cover - simple setter
        self.closed = True
        if self._closed_log is not None:
            self._closed_log.append(self.headers["User-Agent"])


@pytest.fixture(autouse=True)
def reset_pubchem_session() -> None:
    """Restore the shared PubChem session after each test."""

    with pc._SESSION_LOCK:  # type: ignore[attr-defined]
        original_cfg = pc._SESSION_CFG  # type: ignore[attr-defined]
        original_signature = pc._SESSION_SIGNATURE  # type: ignore[attr-defined]
        original_session = pc._session  # type: ignore[attr-defined]
    try:
        yield
    finally:
        with pc._SESSION_LOCK:  # type: ignore[attr-defined]
            current_session = pc._session  # type: ignore[attr-defined]
            pc._SESSION_CFG = original_cfg  # type: ignore[attr-defined]
            pc._SESSION_SIGNATURE = original_signature  # type: ignore[attr-defined]
            pc._session = session_with_retry(*original_cfg)  # type: ignore[attr-defined]
        if current_session is not None:
            current_session.close()


def _configure_session(
    user_agent: str = "pubchem-tests/1.0 (mailto:tests@example.org)",
) -> ApiCfg:
    """Initialise the PubChem session for tests and return the API config."""

    api_cfg = ApiCfg(user_agent=user_agent)
    retry_cfg = RetryCfg()
    pc.init_session(api_cfg, retry_cfg)
    pc.get_session(api_cfg)
    return api_cfg


def test_make_request_serves_cached_results_across_threads(monkeypatch) -> None:
    """Concurrent calls should reuse cached responses safely."""

    url = "https://example.org/data"
    cfg = pl.PubChemCfg(retries=1, delay=0)
    payload = {"ok": True}

    api_cfg = _configure_session()
    session = pc.get_session(api_cfg)

    monkeypatch.setattr(pc, "get_limiter", lambda *args, **kwargs: NoopLimiter())
    monkeypatch.setattr(pc, "sleep", lambda *_: None)
    monkeypatch.setattr(
        session,
        "get",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("cache not used")),
    )

    with pc._CACHE_LOCK:
        pc._CACHE = TTLCache(maxsize=1024, ttl=cfg.cache_ttl)
        pc._CACHE[url] = pc._CacheEntry(payload=payload, outcome="hit")  # type: ignore[attr-defined]

    results: list[dict[str, object]] = []
    start = threading.Event()

    def worker() -> None:
        start.wait()
        result = pc.make_request(url, cfg)
        assert result is payload
        results.append(result)

    threads = [threading.Thread(target=worker) for _ in range(5)]
    for t in threads:
        t.start()
    start.set()
    for t in threads:
        t.join()

    assert len(results) == 5
    with pc._CACHE_LOCK:
        pc._CACHE = None


def test_make_request_cache_reinitialisation_threadsafe(monkeypatch) -> None:
    """Reinitialising the cache under concurrency should be safe."""

    url = "https://example.org/resource"
    payload = {"value": 1}

    calls: list[int] = []
    call_lock = threading.Lock()

    def fake_get(url: str, timeout: tuple[int, int]) -> DummyResponse:
        with call_lock:
            calls.append(threading.get_ident())
        return DummyResponse(payload)

    api_cfg = _configure_session()
    session = pc.get_session(api_cfg)

    monkeypatch.setattr(session, "get", fake_get)
    monkeypatch.setattr(pc, "get_limiter", lambda *args, **kwargs: NoopLimiter())
    monkeypatch.setattr(pc, "sleep", lambda *_: None)

    with pc._CACHE_LOCK:
        pc._CACHE = None

    initial_cfg = pl.PubChemCfg(retries=1, delay=0, cache_ttl=1)
    updated_cfg = pl.PubChemCfg(retries=1, delay=0, cache_ttl=2)

    # Prime the cache with the initial configuration.
    assert pc.make_request(url, initial_cfg) == payload

    def worker(cfg: pl.PubChemCfg) -> dict[str, object] | None:
        return pc.make_request(url, cfg)

    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(worker, cfg) for cfg in (initial_cfg, updated_cfg)]
        results = [f.result() for f in futures]

    assert all(result == payload for result in results)
    assert len(calls) >= 1
    with pc._CACHE_LOCK:
        pc._CACHE = None


def test_make_request_closes_response(monkeypatch) -> None:
    """HTTP responses should be closed after use."""

    url = "https://example.org/close-check"
    payload = {"ok": True}
    cfg = pl.PubChemCfg(retries=0, delay=0)
    close_log: list[object] = []

    api_cfg = _configure_session()
    session = pc.get_session(api_cfg)

    monkeypatch.setattr(pc, "get_limiter", lambda *args, **kwargs: NoopLimiter())
    monkeypatch.setattr(pc, "sleep", lambda *_: None)

    def fake_get(url: str, timeout: tuple[int, int]) -> DummyResponse:
        return DummyResponse(payload, close_log=close_log)

    monkeypatch.setattr(session, "get", fake_get)

    with pc._CACHE_LOCK:
        pc._CACHE = None

    result = pc.make_request(url, cfg)

    assert result == payload
    assert close_log == [True]


def test_get_session_initialises_once_under_concurrency(monkeypatch) -> None:
    """Concurrent access should only create the session a single time."""

    api_cfg = ApiCfg(user_agent="pubchem-tests/2.0 (mailto:tests@example.org)")
    retry_cfg = RetryCfg()
    created: list[DummySession] = []

    def fake_session_with_retry(api: ApiCfg, retry: RetryCfg) -> DummySession:
        session = DummySession(api.user_agent)
        created.append(session)
        return session

    monkeypatch.setattr(pc, "session_with_retry", fake_session_with_retry)
    pc.init_session(api_cfg, retry_cfg)

    with ThreadPoolExecutor(max_workers=8) as pool:
        sessions = list(pool.map(lambda _: pc.get_session(api_cfg), range(8)))

    assert len(created) == 1
    assert len({id(session) for session in sessions}) == 1


def test_init_session_closes_previous_session(monkeypatch) -> None:
    """Reinitialising the session should close the previous instance exactly once."""

    retry_cfg = RetryCfg()
    first_api = ApiCfg(user_agent="pubchem-tests/3.0 (mailto:tests@example.org)")
    second_api = ApiCfg(user_agent="pubchem-tests/4.0 (mailto:tests@example.org)")
    closed: list[str] = []

    def fake_session_with_retry(api: ApiCfg, retry: RetryCfg) -> DummySession:
        return DummySession(api.user_agent, closed)

    monkeypatch.setattr(pc, "session_with_retry", fake_session_with_retry)

    pc.init_session(first_api, retry_cfg)
    first_session = pc.get_session(first_api)
    assert isinstance(first_session, DummySession)
    assert not closed

    pc.init_session(second_api, retry_cfg)
    second_session = pc.get_session(second_api)

    assert isinstance(second_session, DummySession)
    assert closed == [first_api.user_agent]
    assert first_session.closed is True
    assert second_session.closed is False
