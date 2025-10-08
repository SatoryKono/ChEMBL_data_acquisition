from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import pytest

from library.clients import pubchem
from library.config import PubChemCfg


@dataclass
class _DummyResponse:
    status_code: int
    headers: dict[str, Any]
    payload: dict[str, Any] | None = None

    def __enter__(self) -> "_DummyResponse":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: object | None,
    ) -> bool:
        return False

    def json(self) -> dict[str, Any]:
        return self.payload or {}

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise RuntimeError(f"unexpected status {self.status_code}")


class _DummySession:
    def __init__(self, response_factory: Callable[[], _DummyResponse]) -> None:
        self._response_factory = response_factory
        self.calls: list[tuple[str, str, dict[str, Any]]] = []

    def request(
        self, method: str, url: str, **kwargs: Any
    ) -> _DummyResponse:
        self.calls.append((method, url, kwargs))
        return self._response_factory()


class _DummyLimiter:
    def __init__(self) -> None:
        self.acquires: int = 0

    def acquire(self) -> None:
        self.acquires += 1


def test_make_request__caches_server_error_results(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()
    cfg.retries = 1

    url = (
        f"{cfg.base.rstrip('/')}/compound/cid/64972/property/"
        "MolecularFormula,IUPACName,IsomericSMILES,CanonicalSMILES,InChI,InChIKey/JSON"
    )

    sleep_calls: list[float] = []
    monkeypatch.setattr(pubchem, "sleep", lambda seconds: sleep_calls.append(float(seconds)))
    limiter = _DummyLimiter()
    monkeypatch.setattr(pubchem, "get_limiter", lambda *_args, **_kwargs: limiter)

    def _response() -> _DummyResponse:
        return _DummyResponse(503, {"Retry-After": "30"})

    session = _DummySession(_response)
    monkeypatch.setattr(pubchem, "get_session", lambda *_args, **_kwargs: session)
    monkeypatch.setattr(pubchem, "_CACHE", None)

    first_result = pubchem.make_request(url, cfg)

    assert first_result is None
    assert session.calls == [
        (
            "GET",
            url,
            {"timeout": (cfg.timeout_connect, cfg.timeout_read)},
        )
        for _ in range(cfg.retries + 1)
    ]
    assert sleep_calls == [30.0]
    assert limiter.acquires == cfg.retries + 1

    cache = pubchem._ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
    entry = cache.get(pubchem._build_cache_key("GET", url))
    assert entry is not None
    assert entry.outcome == "server_error"
    assert entry.details == {"reason": "server_error", "status": 503, "retry_after": 30.0}

    second_result = pubchem.make_request(url, cfg)

    assert second_result is None
    assert len(session.calls) == cfg.retries + 1
    assert limiter.acquires == cfg.retries + 1
    assert sleep_calls == [30.0]


def test_make_request__invalid_identifier_cached(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = PubChemCfg()

    url = (
        f"{cfg.base.rstrip('/')}/compound/inchi/"
        "InChI%3D1S%2FC7H8N4O2%2Fc8-6-4-2-1-3-5(6)7(9)11-10%2Fh1-4H%2C8H2%2C(,9,10,11)"
        "/cids/JSON"
    )

    limiter = _DummyLimiter()
    monkeypatch.setattr(pubchem, "get_limiter", lambda *_args, **_kwargs: limiter)

    session = _DummySession(lambda: _DummyResponse(400, {}))
    monkeypatch.setattr(pubchem, "get_session", lambda *_args, **_kwargs: session)
    monkeypatch.setattr(pubchem, "_CACHE", None)

    result = pubchem.make_request(url, cfg)

    assert result is None
    assert session.calls == [
        (
            "GET",
            url,
            {"timeout": (cfg.timeout_connect, cfg.timeout_read)},
        )
    ]
    assert limiter.acquires == 1

    cache = pubchem._ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
    entry = cache.get(pubchem._build_cache_key("GET", url))
    assert entry is not None
    assert entry.outcome == "invalid_identifier"
    assert entry.details == {"reason": "invalid_identifier", "status": 400}

    second = pubchem.make_request(url, cfg)

    assert second is None
    assert session.calls == [
        (
            "GET",
            url,
            {"timeout": (cfg.timeout_connect, cfg.timeout_read)},
        )
    ]
    assert limiter.acquires == 1


def test_get_cid_from_smiles__uses_post_for_stereochemistry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()
    smiles = r"C/C=C\C"

    limiter = _DummyLimiter()
    monkeypatch.setattr(pubchem, "get_limiter", lambda *_args, **_kwargs: limiter)

    response_payload = {
        "IdentifierList": {
            "CID": ["123", "123", "456"],
        }
    }

    session = _DummySession(lambda: _DummyResponse(200, {}, response_payload))
    monkeypatch.setattr(pubchem, "get_session", lambda *_args, **_kwargs: session)
    monkeypatch.setattr(pubchem, "_CACHE", None)

    result = pubchem.get_cid_from_smiles(smiles, cfg)

    assert result == "123|456"
    assert session.calls == [
        (
            "POST",
            f"{cfg.base.rstrip('/')}/compound/smiles/cids/JSON",
            {
                "timeout": (cfg.timeout_connect, cfg.timeout_read),
                "data": {"smiles": smiles},
            },
        )
    ]
    assert limiter.acquires == 1

    cache = pubchem._ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
    cache_key = pubchem._build_cache_key(
        "POST",
        f"{cfg.base.rstrip('/')}/compound/smiles/cids/JSON",
        {"smiles": smiles},
    )
    entry = cache.get(cache_key)
    assert entry is not None
    assert entry.is_hit
