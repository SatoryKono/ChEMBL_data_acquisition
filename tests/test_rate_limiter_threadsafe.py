"""Thread-safety tests for :mod:`library.rate_limiter`."""

from __future__ import annotations

import threading

from library import rate_limiter as rl


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
