"""Thread-safety tests for :mod:`library.rate_limiter`."""

from __future__ import annotations

import threading

import pytest

from library import rate_limiter as rl


class ThreadSafeFakeTime:
    """Thread-safe monotonic clock for testing."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._now = 0.0
        self.sleeps: list[float] = []

    def monotonic(self) -> float:
        with self._lock:
            return self._now

    def sleep(self, delay: float) -> None:
        with self._lock:
            self.sleeps.append(delay)
            self._now += delay


def test_get_limiter_thread_safe() -> None:
    """Concurrent calls to :func:`get_limiter` should return the same instance."""
    results: list[rl.RateLimiter] = []

    with rl._limiters_lock:
        rl._limiters.clear()

    def target() -> None:
        results.append(rl.get_limiter("shared", rps=1))

    threads = [threading.Thread(target=target) for _ in range(20)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert len({id(limiter) for limiter in results}) == 1

    with rl._limiters_lock:
        assert len(rl._limiters) == 1
        rl._limiters.clear()


def test_acquire_releases_lock_before_sleep(monkeypatch) -> None:
    """Ensure :meth:`RateLimiter.acquire` sleeps without holding the lock."""

    fake_time = ThreadSafeFakeTime()
    monkeypatch.setattr(rl, "time", fake_time)

    limiter = rl.RateLimiter(rps=1, burst=1)

    def fake_sleep(delay: float) -> None:
        assert not limiter._lock.locked()
        fake_time.sleep(delay)

    monkeypatch.setattr(rl, "sleep", fake_sleep)

    barrier = threading.Barrier(2)
    completed: list[None] = []

    def worker() -> None:
        barrier.wait()
        limiter.acquire()
        completed.append(None)

    threads = [threading.Thread(target=worker) for _ in range(2)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    assert len(completed) == 2
    assert fake_time.sleeps == pytest.approx([1.0])
