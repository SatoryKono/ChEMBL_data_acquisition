"""Tests for :mod:`library.common.rate_limiter`."""

from __future__ import annotations

import threading
import time
from typing import Callable

import pytest

from library.common import rate_limiter


class _FakeClock:
    """Deterministic clock used to fake :func:`time.monotonic` and sleeps."""

    def __init__(self, *, factor: float = 1.0) -> None:
        self._factor = factor
        self._now = 0.0
        self._lock = threading.Lock()
        self.calls: list[float] = []

    def monotonic(self) -> float:
        with self._lock:
            return self._now

    def sleep(self, delay: float) -> None:
        with self._lock:
            self.calls.append(delay)
            self._now += delay * self._factor
        # Yield control to allow other threads to progress in tests using
        # real threading primitives.
        time.sleep(0)

    def reset(self) -> None:
        with self._lock:
            self._now = 0.0
            self.calls.clear()


def _install_clock(monkeypatch: pytest.MonkeyPatch, clock: _FakeClock) -> None:
    """Patch the rate limiter to use *clock* for time keeping."""

    monkeypatch.setattr(rate_limiter.time, "monotonic", clock.monotonic)
    monkeypatch.setattr(rate_limiter, "sleep", clock.sleep)


@pytest.mark.unit
def test_rate_limiter__adaptive_wait_increases_to_cap(monkeypatch: pytest.MonkeyPatch) -> None:
    """The adaptive wait grows until capped when progress is slow."""

    clock = _FakeClock(factor=0.01)
    _install_clock(monkeypatch, clock)

    limiter = rate_limiter.RateLimiter(rps=1.0, burst=1)
    limiter.acquire()  # consume the initial burst token

    clock.reset()
    limiter.acquire()

    assert len(clock.calls) > 3
    assert clock.calls == sorted(clock.calls)
    assert clock.calls[0] == pytest.approx(1.0)
    assert clock.calls[-1] == pytest.approx(rate_limiter._MAX_WAIT_SECONDS)


@pytest.mark.unit
def test_rate_limiter__parallel_threads_make_progress(monkeypatch: pytest.MonkeyPatch) -> None:
    """Concurrent callers all acquire a token without stalling indefinitely."""

    clock = _FakeClock()
    _install_clock(monkeypatch, clock)

    limiter = rate_limiter.RateLimiter(rps=2.0, burst=1)
    acquire_times: list[float] = []
    acquire_lock = threading.Lock()

    def _worker(record: Callable[[float], None]) -> None:
        limiter.acquire()
        record(clock.monotonic())

    def _record(timestamp: float) -> None:
        with acquire_lock:
            acquire_times.append(timestamp)

    threads = [threading.Thread(target=_worker, args=(_record,)) for _ in range(3)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=1)

    assert all(not thread.is_alive() for thread in threads)
    assert len(acquire_times) == 3
    assert acquire_times == sorted(acquire_times)
    assert acquire_times[0] == pytest.approx(0.0)
    assert acquire_times[-1] == pytest.approx(sum(clock.calls))
