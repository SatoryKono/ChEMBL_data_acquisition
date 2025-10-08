"""Simple rate limiting utilities.

This module provides a token-bucket :class:`RateLimiter` and convenience
helpers used across the project to throttle HTTP requests.  Network-facing
code should use :func:`get_limiter` to retrieve a shared limiter instance and
call :meth:`RateLimiter.acquire` before performing a request.  All sleeping is
routed through :func:`sleep` so tests can monkeypatch it easily.
"""

from __future__ import annotations

import math
import threading
import time

from cachetools import TTLCache


GLOBAL_LIMITER_NAME = "system_global"


_MIN_WAIT_SECONDS = 1e-3
_MAX_WAIT_SECONDS = 5.0


class RateLimiter:
    """Token bucket rate limiter.

    Parameters
    ----------
    rps:
        Allowed requests per second.  A value ``<= 0`` disables throttling.
    burst:
        Maximum burst size.  Defaults to ``ceil(rps)``.
    """

    def __init__(self, rps: float, burst: int | None = None) -> None:
        self.rps = rps
        self.burst = burst if burst is not None else max(1, int(rps))
        self._tokens = float(self.burst)
        self._updated = time.monotonic()
        self._lock = threading.Lock()

    def acquire(self) -> None:
        """Block until a token is available.

        The wait strategy uses an adaptive back-off that doubles the minimum
        sleep duration on every retry until a reasonable upper bound is
        reached.  This approach keeps the loop responsive for highly
        contended call sites (including concurrent threads) while guaranteeing
        forward progress even when the system scheduler undersleeps.
        """
        if self.rps <= 0:
            return

        base_wait = max(1.0 / self.rps, _MIN_WAIT_SECONDS)
        max_wait = max(base_wait, _MAX_WAIT_SECONDS)
        adaptive_wait = base_wait
        wait = 0.0

        while True:
            if wait > 0:
                sleep(wait)

            with self._lock:
                now = time.monotonic()
                elapsed = now - self._updated
                self._tokens = min(float(self.burst), self._tokens + elapsed * self.rps)
                if self._tokens >= 1 or math.isclose(
                    self._tokens, 1.0, rel_tol=0.0, abs_tol=1e-9
                ):
                    self._tokens = max(0.0, self._tokens - 1)
                    self._updated = now
                    return

                required = max(0.0, (1 - self._tokens) / self.rps)
                wait = min(max(required, adaptive_wait), max_wait)

            adaptive_wait = min(adaptive_wait * 2, max_wait)


_limiters: TTLCache[str, RateLimiter] = TTLCache(maxsize=128, ttl=600)

_limiters_lock = threading.Lock()


def configure_limiter_cache(maxsize: int, ttl: int) -> None:
    """Configure the cache storing :class:`RateLimiter` instances.

    Parameters
    ----------
    maxsize:
        Maximum number of cached limiters.
    ttl:
        Time-to-live for cache entries in seconds.
    """

    global _limiters
    with _limiters_lock:
        _limiters = TTLCache(maxsize=maxsize, ttl=ttl)


def get_limiter(name: str, rps: float, burst: int | None = None) -> RateLimiter:
    """Return a shared :class:`RateLimiter` identified by ``name``.

    Parameters
    ----------
    name:
        Identifier for the limiter.  Subsequent calls with the same name return
        the same instance.
    rps:
        Allowed requests per second.  A value ``<= 0`` disables throttling.
    burst:
        Maximum burst size.  Defaults to ``ceil(rps)``.
    """
    with _limiters_lock:
        limiter: RateLimiter | None = _limiters.get(name)

        if (
            limiter is None
            or limiter.rps != rps
            or (burst is not None and limiter.burst != burst)
        ):
            limiter = RateLimiter(rps, burst)
            _limiters[name] = limiter
        return limiter


def get_global_limiter(
    rps: float | None, burst: int | None = None
) -> RateLimiter | None:
    """Return the shared system-wide :class:`RateLimiter` if enabled.

    Parameters
    ----------
    rps:
        Allowed requests per second for the entire pipeline.  A value ``<= 0``
        disables the limiter.
    burst:
        Maximum burst size for the global limiter.  Non-positive values are
        treated as ``None``.
    """

    if rps is None or rps <= 0:
        return None
    burst_value = burst if burst is not None and burst > 0 else None
    return get_limiter(GLOBAL_LIMITER_NAME, rps, burst_value)


def sleep(delay: float) -> None:
    """Sleep for ``delay`` seconds.

    Parameters
    ----------
    delay:
        Number of seconds to pause execution.
    """
    time.sleep(delay)
