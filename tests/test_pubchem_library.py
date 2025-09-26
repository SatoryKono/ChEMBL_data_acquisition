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


@responses.activate
def test_resolve_pubchem_record_falls_back_to_name() -> None:
    """The resolver should follow ``resolve_order`` when lookups fail."""

    cfg = pl.PubChemCfg(
        resolve_order=("smiles", "pref_name"),
        retries=1,
        delay=0,
        cache_ttl=0,
    )

    smiles_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/C/"
        "property/MolecularFormula,IUPACName,IsomericSMILES,CanonicalSMILES,InChI,InChIKey/JSON"
    )
    name_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Aspirin/"
        "property/MolecularFormula,IUPACName,IsomericSMILES,CanonicalSMILES,InChI,InChIKey/JSON"
    )

    responses.add(responses.GET, smiles_url, status=404)
    responses.add(
        responses.GET,
        name_url,
        json={
            "PropertyTable": {
                "Properties": [
                    {
                        "CID": 2244,
                        "IUPACName": "Aspirin",
                        "MolecularFormula": "C9H8O4",
                        "IsomericSMILES": "O=C(O)C1=CC=CC=C1OC(=O)C",
                        "CanonicalSMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
                        "InChI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
                        "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                    }
                ]
            }
        },
        status=200,
    )

    record = pl.resolve_pubchem_record(
        {"canonical_smiles": "C", "pref_name": "Aspirin"},
        cfg,
        cache={},
    )

    assert record["pubchem_cid"] == "2244"
    assert record["pubchem_iupac_name"] == "Aspirin"
    assert len(responses.calls) == 2
    assert responses.calls[0].request.url == smiles_url
    assert responses.calls[1].request.url == name_url


@responses.activate
def test_resolve_pubchem_record_backoff_on_status() -> None:
    """HTTP 429 and 5xx responses should trigger exponential backoff."""

    cfg = pl.PubChemCfg(
        resolve_order=("smiles",),
        retries=3,
        delay=0,
        cache_ttl=0,
        backoff_initial_seconds=0.1,
        timeout_seconds=5,
    )

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/C/"
        "property/MolecularFormula,IUPACName,IsomericSMILES,CanonicalSMILES,InChI,InChIKey/JSON"
    )

    responses.add(responses.GET, url, status=429)
    responses.add(responses.GET, url, status=503)
    responses.add(
        responses.GET,
        url,
        json={
            "PropertyTable": {
                "Properties": [
                    {
                        "CID": 5793,
                        "IUPACName": "Acetic acid",
                        "MolecularFormula": "C2H4O2",
                        "IsomericSMILES": "CC(=O)O",
                        "CanonicalSMILES": "CC(=O)O",
                        "InChI": "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)",
                        "InChIKey": "QTBSBXVTEAMEQO-UHFFFAOYSA-N",
                    }
                ]
            }
        },
        status=200,
    )

    sleeps: list[float] = []

    def fake_sleep(delay: float) -> None:
        sleeps.append(delay)

    record = pl.resolve_pubchem_record(
        {"canonical_smiles": "C"},
        cfg,
        cache={},
        sleeper=fake_sleep,
    )

    assert record["pubchem_cid"] == "5793"
    assert len(responses.calls) == 3
    assert sleeps == [0.1, 0.2]
