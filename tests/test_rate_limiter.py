import time

from library.rate_limiter import RateLimiter


def test_rate_limiter_enforces_rps() -> None:
    """Ensure the limiter respects the configured requests per second."""
    limiter = RateLimiter(rps=5, burst=1)
    n_calls = 5
    start = time.monotonic()
    for _ in range(n_calls):
        limiter.acquire()
    elapsed = time.monotonic() - start
    assert elapsed >= (n_calls - 1) / 5
