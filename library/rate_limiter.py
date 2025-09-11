from __future__ import annotations

import threading
import time

from .config import RateCfg, load_config


class RateLimiter:
    """Token-bucket rate limiter.

    Parameters
    ----------
    rps:
        Target requests per second.
    burst:
        Maximum burst size allowed.
    """

    def __init__(self, rps: float, burst: int) -> None:
        self._rps = float(rps)
        self._capacity = float(max(1, burst))
        self._tokens = self._capacity
        self._updated = time.monotonic()
        self._lock = threading.Lock()

    def wait(self, extra_delay: float = 0.0) -> None:
        """Block until a token is available and sleep ``extra_delay`` seconds.

        Parameters
        ----------
        extra_delay:
            Additional time in seconds to wait after acquiring a token.
        """
        if self._rps <= 0:
            if extra_delay > 0:
                time.sleep(extra_delay)
            return
        while True:
            with self._lock:
                now = time.monotonic()
                elapsed = now - self._updated
                self._tokens = min(self._capacity, self._tokens + elapsed * self._rps)
                if self._tokens >= 1:
                    self._tokens -= 1
                    self._updated = now
                    break
                wait_time = (1 - self._tokens) / self._rps
                self._updated = now
            time.sleep(wait_time)
        if extra_delay > 0:
            time.sleep(extra_delay)


try:
    _cfg = load_config()
    rate_limiter = RateLimiter(_cfg.rate.global_rps, _cfg.rate.global_burst)
except Exception:
    _defaults = RateCfg()
    rate_limiter = RateLimiter(_defaults.global_rps, _defaults.global_burst)


__all__ = ["RateLimiter", "rate_limiter"]
