from library import rate_limiter as rl


class FakeTime:
    def __init__(self) -> None:
        self.now = 0.0
        self.sleeps: list[float] = []

    def monotonic(self) -> float:
        return self.now

    def sleep(self, delay: float) -> None:
        self.sleeps.append(delay)
        self.now += delay


def test_rate_limiter_enforces_rps(monkeypatch) -> None:
    """Ensure the limiter respects the configured requests per second."""
    fake_time = FakeTime()
    monkeypatch.setattr(rl, "time", fake_time)
    monkeypatch.setattr(rl, "sleep", fake_time.sleep)
    limiter = rl.RateLimiter(rps=5, burst=1)
    for _ in range(5):
        limiter.acquire()
    assert fake_time.sleeps == [0.2, 0.2, 0.2, 0.2]


def test_rate_limiter_combined_throttling(monkeypatch) -> None:
    """Global and service limiters should compose into a combined delay."""

    fake_time = FakeTime()
    monkeypatch.setattr(rl, "time", fake_time)
    monkeypatch.setattr(rl, "sleep", fake_time.sleep)

    documents = rl.RateLimiter(rps=10, burst=1)
    service = rl.RateLimiter(rps=4, burst=1)
    request = rl.RateLimiter(rps=2, burst=1)

    for _ in range(3):
        documents.acquire()
        service.acquire()
        request.acquire()

    assert fake_time.sleeps == [0.1, 0.25, 0.5, 0.5]


def test_get_global_limiter_disabled_returns_none() -> None:
    """Zero rates should avoid creating a cached global limiter."""

    with rl._limiters_lock:
        rl._limiters.clear()

    limiter = rl.get_global_limiter(0, 0)

    assert limiter is None
    with rl._limiters_lock:
        assert rl.GLOBAL_LIMITER_NAME not in rl._limiters
