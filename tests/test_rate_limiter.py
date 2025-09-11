import time

from library.rate_limiter import RateLimiter


def test_rate_limiter_enforces_rps() -> None:
    """Rate limiter should enforce the configured requests per second."""
    limiter = RateLimiter(2, 1)
    start = time.perf_counter()
    for _ in range(3):
        limiter.wait()
    elapsed = time.perf_counter() - start
    assert elapsed >= 1.0
