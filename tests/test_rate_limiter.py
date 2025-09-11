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
