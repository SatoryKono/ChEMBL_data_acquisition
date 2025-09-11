"""Tests for limiter cache size in :mod:`library.rate_limiter`."""

from __future__ import annotations

import library.rate_limiter as rl


def test_limiter_cache_respects_maxsize() -> None:
    """After many limiter requests, cache size should not exceed ``maxsize``."""
    rl.configure_limiter_cache(maxsize=100, ttl=60)
    for i in range(10_000):
        rl.get_limiter(f"name{i}", rps=1)
    with rl._limiters_lock:
        assert len(rl._limiters) <= 100
    rl.configure_limiter_cache(maxsize=128, ttl=600)
