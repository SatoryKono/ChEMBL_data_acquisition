"""Tests for :mod:`library.pubchem_library`."""

from __future__ import annotations

import time

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
    pl._CACHE = None

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
    pl._CACHE = None

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
    pl._CACHE = None

    cfg = pl.PubChemCfg(retries=2, delay=1)
    pl.make_request("https://example.org", cfg)

    assert sleeps == [1]


@responses.activate
def test_get_properties_returns_none_for_missing() -> None:
    """Missing PubChem fields should be returned as ``None``."""

    cfg = pl.PubChemCfg(delay=0)
    cid = "123"
    url = (
        f"{cfg.base.rstrip('/')}/compound/cid/{cid}/property/"
        "MolecularFormula,IUPACName,IsomericSMILES,CanonicalSMILES,InChI,InChIKey/JSON"
    )
    responses.add(responses.GET, url, json={"PropertyTable": {"Properties": [{}]}}, status=200)

    props = pl.get_properties(cid, cfg)

    assert props.IUPACName is None
    assert props.MolecularFormula is None
    assert props.iSMILES is None
    assert props.cSMILES is None
    assert props.InChI is None
    assert props.InChIKey is None


def test_cache_entry_expires(monkeypatch: pytest.MonkeyPatch) -> None:
    """Cache entries should be evicted after the configured TTL."""

    class Resp:
        status_code = 200

        def json(self) -> dict[str, object]:
            return {}

        def raise_for_status(self) -> None:  # pragma: no cover - no error
            return None

    calls = {"n": 0}

    def capture(url: str, timeout: tuple[int, int]) -> Resp:
        calls["n"] += 1
        return Resp()

    monkeypatch.setattr(pl._session, "get", capture)
    monkeypatch.setattr(
        pl,
        "get_limiter",
        lambda *a, **k: type("L", (), {"acquire": lambda self: None})(),
    )
    pl._CACHE = None

    cfg = pl.PubChemCfg(cache_ttl=1, delay=0, retries=1)
    url = "https://example.org"

    pl.make_request(url, cfg)
    pl.make_request(url, cfg)
    assert calls["n"] == 1  # cached

    time.sleep(1.1)
    pl.make_request(url, cfg)
    assert calls["n"] == 2  # cache expired


def test_resolve_pubchem_record_follows_order(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = pl.PubChemCfg(delay=0, retries=1, backoff_initial_seconds=0.1)

    calls: list[tuple[str, bool]] = []

    def fake_make_request(
        url: str,
        cfg: pl.PubChemCfg,
        *,
        deadline: float | None = None,
        return_status: bool = False,
    ) -> dict[str, object] | tuple[int | None, dict[str, object] | None] | None:
        calls.append((url, return_status))
        if "smiles" in url:
            if return_status:
                return 404, None
            return None
        if "inchikey" in url:
            data = {"IdentifierList": {"CID": ["321"]}}
            if return_status:
                return 200, data
            return data
        if "property" in url:
            return {"PropertyTable": {"Properties": [{"IUPACName": "foo"}]}}
        raise AssertionError(f"unexpected URL {url}")

    monkeypatch.setattr(pl, "make_request", fake_make_request)
    monkeypatch.setattr(pl, "_CACHE", None)

    cache: dict[tuple[str, str], tuple[bool, dict[str, str | None] | None]] = {}
    identifiers = {
        "cid": None,
        "smiles": "C",
        "inchikey": "AAA",
        "inchi": None,
        "pref_name": None,
    }

    record = pl.resolve_pubchem_record(identifiers, cfg, cache=cache)

    assert record["pubchem_cid"] == "321"
    assert record["pubchem_iupac_name"] == "foo"
    assert calls[0][0].endswith("/smiles/C/cids/JSON")
    assert calls[1][0].endswith("/inchikey/AAA/cids/JSON")
    assert cache[("smiles", "C")] == (False, None)
    found, stored = cache[("inchikey", "AAA")]
    assert found
    assert stored is not None

    call_count = len(calls)
    second = pl.resolve_pubchem_record(identifiers, cfg, cache=cache)
    assert second == record
    assert len(calls) == call_count


def test_make_request_backoff_on_429(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = pl.PubChemCfg(delay=0, retries=2, backoff_initial_seconds=0.25)

    class StubResponse:
        def __init__(self, status_code: int, payload: dict[str, object] | None = None) -> None:
            self.status_code = status_code
            self._payload = payload or {}

        def json(self) -> dict[str, object]:
            return self._payload

        def raise_for_status(self) -> None:
            if self.status_code >= 400:
                raise requests.HTTPError(response=self)

    responses_iter = iter(
        [
            StubResponse(429, {}),
            StubResponse(200, {"IdentifierList": {"CID": [1]}}),
        ]
    )

    def fake_get(url: str, timeout: tuple[int, int]) -> StubResponse:
        return next(responses_iter)

    sleeps: list[float] = []

    monkeypatch.setattr(pl, "sleep", lambda delay: sleeps.append(delay))
    monkeypatch.setattr(
        pl,
        "get_limiter",
        lambda *a, **k: type("L", (), {"acquire": lambda self: None})(),
    )
    monkeypatch.setattr(pl, "_CACHE", None)
    monkeypatch.setattr(pl._session, "get", fake_get)

    status, data = pl.make_request("https://example.org", cfg, return_status=True)

    assert status == 200
    assert data == {"IdentifierList": {"CID": [1]}}
    assert sleeps == [0.25]
