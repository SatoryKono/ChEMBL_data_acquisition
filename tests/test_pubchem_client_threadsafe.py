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

    def __init__(
        self,
        payload: dict[str, object] | None,
        close_log: list[object] | None = None,
        *,
        status_code: int = 200,
        headers: dict[str, str] | None = None,
    ) -> None:
        self._payload = payload
        self._close_log = close_log
        self.status_code = status_code
        self.headers = headers or {}
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
        if self._payload is None:
            raise ValueError("response payload not set")
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


class FakeClock:
    """Deterministic clock used to track sleeps without real delays."""

    def __init__(self) -> None:
        self.current = 0.0
        self.sleeps: list[float] = []

    def monotonic(self) -> float:
        return self.current

    def sleep(self, duration: float) -> None:
        self.sleeps.append(duration)
        self.current += duration


@pytest.fixture(autouse=True)
def reset_pubchem_session() -> None:
    """Restore the shared PubChem session after each test."""

    with pc._SESSION_LOCK:  # type: ignore[attr-defined]
        original_cfg = pc._SESSION_CFG  # type: ignore[attr-defined]
        original_signature = pc._SESSION_SIGNATURE  # type: ignore[attr-defined]
        original_initialised = pc._SESSION_INITIALISED  # type: ignore[attr-defined]
        original_session = pc._session  # type: ignore[attr-defined]
    try:
        yield
    finally:
        with pc._SESSION_LOCK:  # type: ignore[attr-defined]
            current_session = pc._session  # type: ignore[attr-defined]
            pc._SESSION_CFG = original_cfg  # type: ignore[attr-defined]
            pc._SESSION_SIGNATURE = original_signature  # type: ignore[attr-defined]
            pc._SESSION_INITIALISED = original_initialised  # type: ignore[attr-defined]
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


def test_get_session_initialises_for_default_user_agent() -> None:
    """Default configuration should initialise lazily when unused."""

    api_cfg = ApiCfg()

    session = pc.get_session(api_cfg)

    assert session.headers.get("User-Agent") == api_cfg.user_agent


def test_get_session_allows_default_after_initialisation() -> None:
    """Calling init_session enables reuse of the bundled default agent."""

    api_cfg = ApiCfg()
    retry_cfg = RetryCfg()
    pc.init_session(api_cfg, retry_cfg)

    session = pc.get_session()

    assert session.headers.get("User-Agent") == api_cfg.user_agent


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


def test_make_request_refreshes_session_for_new_user_agent(monkeypatch) -> None:
    """Requests with different agents should recreate the underlying session."""

    url_one = "https://example.org/first"
    url_two = "https://example.org/second"
    cfg_one = pl.PubChemCfg(
        retries=0,
        delay=0,
        user_agent="pubchem-tests/4.0 (mailto:tests@example.org)",
    )
    cfg_two = pl.PubChemCfg(
        retries=0,
        delay=0,
        user_agent="pubchem-tests/5.0 (mailto:tests@example.org)",
    )

    created_sessions: list[DummySession] = []
    close_log: list[str] = []

    def fake_session_with_retry(api_cfg: ApiCfg, retry_cfg: RetryCfg) -> DummySession:
        session = DummySession(api_cfg.user_agent, closed_log=close_log)

        def fake_get(url: str, timeout: tuple[int, int]) -> DummyResponse:
            return DummyResponse({"user_agent": session.headers["User-Agent"]})

        session.get = fake_get  # type: ignore[attr-defined]
        created_sessions.append(session)
        return session

    monkeypatch.setattr(pc, "session_with_retry", fake_session_with_retry)
    monkeypatch.setattr(pc, "get_limiter", lambda *args, **kwargs: NoopLimiter())
    monkeypatch.setattr(pc, "sleep", lambda *_: None)

    pc.init_session(
        ApiCfg(user_agent="pubchem-tests/init (mailto:tests@example.org)"), RetryCfg()
    )
    with pc._CACHE_LOCK:
        pc._CACHE = None

    first_response = pc.make_request(url_one, cfg_one)
    second_response = pc.make_request(url_two, cfg_two)

    assert created_sessions[0].headers["User-Agent"] == cfg_one.user_agent
    assert created_sessions[1].headers["User-Agent"] == cfg_two.user_agent
    assert close_log == [cfg_one.user_agent]
    assert first_response == {"user_agent": cfg_one.user_agent}
    assert second_response == {"user_agent": cfg_two.user_agent}


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

    initial_cfg = pl.PubChemCfg(retries=1, delay=0, cache_ttl=1).model_copy(
        update={"user_agent": api_cfg.user_agent}
    )
    updated_cfg = pl.PubChemCfg(retries=1, delay=0, cache_ttl=2).model_copy(
        update={"user_agent": api_cfg.user_agent}
    )

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
    close_log: list[object] = []

    api_cfg = _configure_session()
    session = pc.get_session(api_cfg)

    cfg = pl.PubChemCfg(retries=0, delay=0).model_copy(
        update={"user_agent": api_cfg.user_agent}
    )

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


def test_make_request_retry_after_adjusts_delay(monkeypatch) -> None:
    """Retry-After headers should dictate the backoff delay."""

    url = "https://example.org/rate-limit"
    payload = {"ok": True}

    api_cfg = _configure_session()
    session = pc.get_session(api_cfg)

    cfg = pl.PubChemCfg(retries=1, delay=0, backoff_initial_seconds=0.1).model_copy(
        update={"user_agent": api_cfg.user_agent}
    )

    responses = [
        DummyResponse(
            None,
            status_code=429,
            headers={"Retry-After": "2.5"},
        ),
        DummyResponse(payload),
    ]

    def fake_get(url: str, timeout: tuple[int, int]) -> DummyResponse:
        return responses.pop(0)

    clock = FakeClock()

    monkeypatch.setattr(session, "get", fake_get)
    monkeypatch.setattr(pc, "get_limiter", lambda *args, **kwargs: NoopLimiter())
    monkeypatch.setattr(pc, "sleep", clock.sleep)
    monkeypatch.setattr(pc, "monotonic", clock.monotonic)

    with pc._CACHE_LOCK:
        pc._CACHE = None

    result = pc.make_request(url, cfg)

    assert result == payload
    assert responses == []
    assert clock.sleeps == [pytest.approx(2.5)]

    with pc._CACHE_LOCK:
        pc._CACHE = None


def test_make_request_timeout_records_retry_after(monkeypatch) -> None:
    """Timeouts should cache the last failure details including Retry-After."""

    url = "https://example.org/timeout"
    api_cfg = _configure_session()
    session = pc.get_session(api_cfg)

    cfg = pl.PubChemCfg(
        retries=3,
        delay=0,
        backoff_initial_seconds=0,
        timeout_seconds=3,
    ).model_copy(update={"user_agent": api_cfg.user_agent})

    responses = [
        DummyResponse(None, status_code=429, headers={"Retry-After": "2"}),
        DummyResponse(None, status_code=429, headers={"Retry-After": "2"}),
    ]

    def fake_get(url: str, timeout: tuple[int, int]) -> DummyResponse:
        return responses.pop(0)

    clock = FakeClock()

    monkeypatch.setattr(session, "get", fake_get)
    monkeypatch.setattr(pc, "get_limiter", lambda *args, **kwargs: NoopLimiter())
    monkeypatch.setattr(pc, "sleep", clock.sleep)
    monkeypatch.setattr(pc, "monotonic", clock.monotonic)

    with pc._CACHE_LOCK:
        pc._CACHE = None

    result = pc.make_request(url, cfg)

    assert result is None
    assert responses == []
    assert clock.sleeps == [pytest.approx(2.0), pytest.approx(2.0)]

    with pc._CACHE_LOCK:
        cache = pc._CACHE
        assert cache is not None
        entry = cache[url]

    assert entry.outcome == "timeout"
    assert entry.details is not None
    assert entry.details.get("reason") == "rate_limited"
    assert entry.details.get("status") == 429
    assert entry.details.get("retry_after") == pytest.approx(2.0)

    with pc._CACHE_LOCK:
        pc._CACHE = None


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
