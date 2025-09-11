"""Simple rate limiting utilities.

This module provides a token-bucket :class:`RateLimiter` and convenience
helpers used across the project to throttle HTTP requests.  Network-facing
code should use :func:`get_limiter` to retrieve a shared limiter instance and
call :meth:`RateLimiter.acquire` before performing a request.  All sleeping is
routed through :func:`sleep` so tests can monkeypatch it easily.
"""

from __future__ import annotations

import threading
import time
from typing import Dict


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
        """Block until a token is available."""
        if self.rps <= 0:
            return
        with self._lock:
            now = time.monotonic()
            elapsed = now - self._updated
            self._tokens = min(float(self.burst), self._tokens + elapsed * self.rps)
            if self._tokens < 1:
                wait = (1 - self._tokens) / self.rps
                sleep(wait)
                now = time.monotonic()
                elapsed = now - self._updated
                self._tokens = min(float(self.burst), self._tokens + elapsed * self.rps)
            self._tokens -= 1
            self._updated = now


_limiters: Dict[str, RateLimiter] = {}
_limiters_lock = threading.Lock()


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
        limiter = _limiters.get(name)
        if (
            limiter is None
            or limiter.rps != rps
            or (burst is not None and limiter.burst != burst)
        ):
            limiter = RateLimiter(rps, burst)
            _limiters[name] = limiter
        return limiter


def sleep(delay: float) -> None:
    """Sleep for ``delay`` seconds.

    Parameters
    ----------
    delay:
        Number of seconds to pause execution.
    """
    time.sleep(delay)
