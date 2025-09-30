"""Thread-safety tests for :mod:`library.clients.pubchem`."""

from __future__ import annotations

import threading
from concurrent.futures import ThreadPoolExecutor

from cachetools import TTLCache

from library import pubchem_library as pl
from library.clients import pubchem as pc


class DummyResponse:
    """Simple response object returning a constant payload."""

    status_code = 200

    def __init__(self, payload: dict[str, object]) -> None:
        self._payload = payload

    def __enter__(self) -> DummyResponse:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # pragma: no cover - no-op
        return None

    def raise_for_status(self) -> None:  # pragma: no cover - no error
        return None

    def json(self) -> dict[str, object]:
        return self._payload


class NoopLimiter:
    """Limiter stub avoiding waits in tests."""

    def acquire(self) -> None:  # pragma: no cover - simple stub
        return None


def test_make_request_serves_cached_results_across_threads(monkeypatch) -> None:
    """Concurrent calls should reuse cached responses safely."""

    url = "https://example.org/data"
    cfg = pl.PubChemCfg(retries=1, delay=0)
    payload = {"ok": True}

    monkeypatch.setattr(pc, "get_limiter", lambda *args, **kwargs: NoopLimiter())
    monkeypatch.setattr(pc, "sleep", lambda *_: None)
    monkeypatch.setattr(
        pc._session,  # type: ignore[attr-defined]
        "get",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("cache not used")),
    )

    with pc._CACHE_LOCK:
        pc._CACHE = TTLCache(maxsize=1024, ttl=cfg.cache_ttl)
        pc._CACHE[url] = payload

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

    monkeypatch.setattr(pc._session, "get", fake_get)
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
