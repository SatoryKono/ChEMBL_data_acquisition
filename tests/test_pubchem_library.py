"""Tests for :mod:`library.pubchem_library`."""

from __future__ import annotations

import pytest
import requests

responses = pytest.importorskip("responses")


from library import pubchem_library as pl  # noqa: E402
from library import rate_limiter as rl  # noqa: E402


@responses.activate
def test_get_cid_from_smiles_uses_base() -> None:
    """Ensure the configured base URL is used for PubChem requests."""
    cfg = pl.PubChemCfg(base="https://example.org/api", delay=0)

    url = "https://example.org/api/compound/smiles/C/cids/JSON"
    responses.add(
        responses.GET,
        url,
        json={"IdentifierList": {"CID": [1]}},
        status=200,
    )

    cid = pl.get_cid_from_smiles("C", cfg)

    assert responses.calls[0].request.url == url
    assert cid == "1"


@responses.activate
def test_make_request_uses_timeout(monkeypatch: pytest.MonkeyPatch) -> None:
    """`make_request` passes configured timeouts to the session."""
    called: dict[str, tuple[int, int]] = {}

    class Resp:
        status_code = 200

        def json(self) -> dict[str, object]:
            return {}

        def raise_for_status(self) -> None:  # pragma: no cover - no error
            return None

    def capture(url: str, timeout: tuple[int, int]) -> Resp:
        called["timeout"] = timeout
        return Resp()

    monkeypatch.setattr(pl._session, "get", capture)
    monkeypatch.setattr(
        pl,
        "get_limiter",
        lambda *a, **k: type("L", (), {"acquire": lambda self: None})(),
    )
    pl._CACHE.clear()

    cfg = pl.PubChemCfg(timeout_connect=1, timeout_read=2, retries=1, rps=1)

    responses.add(responses.GET, "https://example.org", json={}, status=200)
    pl.make_request("https://example.org", cfg)
    assert called["timeout"] == (1, 2)


def test_make_request_rate_limited(monkeypatch: pytest.MonkeyPatch) -> None:
    """``make_request`` should pause when RPS is exceeded."""

    class FakeTime:
        def __init__(self) -> None:
            self.current = 0.0
            self.sleeps: list[float] = []

        def monotonic(self) -> float:
            return self.current

        def sleep(self, delay: float) -> None:
            self.sleeps.append(delay)
            self.current += delay

    fake_time = FakeTime()
    monkeypatch.setattr(rl, "time", fake_time)
    monkeypatch.setattr(rl, "sleep", fake_time.sleep)
    with rl._limiters_lock:
        rl._limiters.clear()

    class Resp:
        status_code = 200

        def json(self) -> dict[str, object]:
            return {}

        def raise_for_status(self) -> None:  # pragma: no cover - no error
            return None

    monkeypatch.setattr(pl._session, "get", lambda url, timeout: Resp())
    pl._CACHE.clear()

    cfg = pl.PubChemCfg(rps=1, burst=1, retries=1)
    pl.make_request("https://example.org/a", cfg)
    pl.make_request("https://example.org/b", cfg)

    called: dict[str, tuple[int, int]] = {}

    def capture(url: str, timeout: tuple[int, int]) -> Resp:
        called["timeout"] = timeout
        return Resp()

    monkeypatch.setattr(pl._session, "get", capture)
    cfg = pl.PubChemCfg(timeout_connect=1, timeout_read=2, delay=0, retries=1)

    pl.make_request("https://example.org", cfg)

    assert called["timeout"] == (1, 2)

    assert fake_time.sleeps == [1.0]


def test_make_request_waits_between_retries(monkeypatch) -> None:
    """``make_request`` should delay between failed attempts."""

    sleeps: list[float] = []

    def fake_sleep(delay: float) -> None:
        sleeps.append(delay)

    class Resp:
        status_code = 200

        def json(self) -> dict[str, object]:
            return {}

        def raise_for_status(self) -> None:  # pragma: no cover - no error
            return None

    attempts = {"n": 0}

    def fake_get(url: str, timeout: tuple[int, int]) -> Resp:
        attempts["n"] += 1
        if attempts["n"] == 1:
            raise requests.RequestException("boom")
        return Resp()

    monkeypatch.setattr(pl, "sleep", fake_sleep)
    monkeypatch.setattr(pl._session, "get", fake_get)

    class Limiter:
        def acquire(self) -> None:  # pragma: no cover - simple stub
            return None

    monkeypatch.setattr(pl, "get_limiter", lambda *a, **k: Limiter())
    pl._CACHE.clear()

    cfg = pl.PubChemCfg(retries=2, delay=1)
    pl.make_request("https://example.org", cfg)

    assert sleeps == [1]
