from __future__ import annotations

import random
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, Generator


import pytest

from library.clients import pubchem
from library.config import ApiCfg, PubChemCfg, RetryCfg


@pytest.fixture(autouse=True)
def _reset_service_outage(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(pubchem, "_SERVICE_OUTAGE_UNTIL", None)
    monkeypatch.setattr(pubchem, "_SERVICE_OUTAGE_REASON", None)
    monkeypatch.setattr(pubchem, "_SERVICE_OUTAGE_DETAILS", None)


@pytest.fixture(autouse=True)
def _reset_pubchem_session() -> Generator[None, None, None]:
    original_session = pubchem._session
    original_cfg = pubchem._SESSION_CFG
    original_signature = pubchem._SESSION_SIGNATURE
    original_initialised = pubchem._SESSION_INITIALISED

    yield

    restored_session = pubchem._session
    if restored_session is not None and restored_session is not original_session:
        restored_session.close()
    pubchem._session = original_session
    pubchem._SESSION_CFG = original_cfg
    pubchem._SESSION_SIGNATURE = original_signature
    pubchem._SESSION_INITIALISED = original_initialised


def test_init_session__rejects_placeholder_user_agent() -> None:
    api_cfg = ApiCfg.model_construct(user_agent="contact@example.org")
    retry_cfg = RetryCfg()

    with pytest.raises(ValueError, match="User-Agent"):
        pubchem.init_session(api_cfg, retry_cfg)


def test_get_session__rejects_placeholder_user_agent() -> None:
    valid_api = ApiCfg()
    retry_cfg = RetryCfg()

    pubchem.init_session(valid_api, retry_cfg)

    with pytest.raises(ValueError, match="User-Agent"):
        pubchem.get_session(ApiCfg.model_construct(user_agent="contact@example.org"))


@dataclass
class _DummyResponse:
    status_code: int
    headers: dict[str, Any]
    payload: dict[str, Any] | None = None

    def __enter__(self) -> _DummyResponse:
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

    def request(self, method: str, url: str, **kwargs: Any) -> _DummyResponse:
        self.calls.append((method, url, kwargs))
        return self._response_factory()


class _DummyLimiter:
    def __init__(self) -> None:
        self.acquires: int = 0

    def acquire(self) -> None:
        self.acquires += 1


def test_make_request__applies_verify_setting(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = PubChemCfg(verify="/tmp/custom-ca.pem")
    cfg.retries = 0

    payload = {"status": "ok"}
    url = f"{cfg.base.rstrip('/')}/compound/cid/42/property/Foo/JSON"

    limiter = _DummyLimiter()
    monkeypatch.setattr(pubchem, "get_limiter", lambda *_args, **_kwargs: limiter)

    session = _DummySession(lambda: _DummyResponse(200, {}, payload))
    session.verify = True  # type: ignore[attr-defined]
    monkeypatch.setattr(pubchem, "get_session", lambda *_args, **_kwargs: session)
    monkeypatch.setattr(pubchem, "_CACHE", None)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_UNTIL", None)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_DETAILS", None)

    result = pubchem.make_request(url, cfg)

    assert result == payload
    assert session.verify == cfg.verify  # type: ignore[attr-defined]
    assert session.calls == [
        (
            "GET",
            url,
            {"timeout": (cfg.timeout_connect, cfg.timeout_read)},
        )
    ]
    assert limiter.acquires == 1


class _TimeController:
    def __init__(self) -> None:
        self._now = 0.0

    def monotonic(self) -> float:
        return self._now

    def sleep(self, seconds: float) -> None:
        self._now += float(seconds)

    def advance(self, seconds: float) -> None:
        self._now += float(seconds)


def test_make_request__caches_server_error_results(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()
    cfg.retries = 1
    cfg.timeout_seconds = 120

    url = (
        f"{cfg.base.rstrip('/')}/compound/cid/64972/property/"
        "MolecularFormula,IUPACName,IsomericSMILES,CanonicalSMILES,InChI,InChIKey/JSON"
    )

    sleep_calls: list[float] = []
    monkeypatch.setattr(
        pubchem, "sleep", lambda seconds: sleep_calls.append(float(seconds))
    )
    limiter = _DummyLimiter()
    monkeypatch.setattr(pubchem, "get_limiter", lambda *_args, **_kwargs: limiter)

    def _response() -> _DummyResponse:
        return _DummyResponse(503, {"Retry-After": "30"})

    session = _DummySession(_response)
    monkeypatch.setattr(pubchem, "get_session", lambda *_args, **_kwargs: session)
    monkeypatch.setattr(pubchem, "_CACHE", None)

    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_UNTIL", None)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_DETAILS", None)

    first_result = pubchem.make_request(url, cfg)

    assert first_result is None
    assert session.calls == [
        (
            "GET",
            url,
            {"timeout": (cfg.timeout_connect, cfg.timeout_read)},
        )
    ]
    assert sleep_calls == []
    assert limiter.acquires == 1
    outcome, details = pubchem.last_request_outcome()
    assert outcome == "server_error"
    assert details is not None
    assert details.get("reason") == "server_error"
    assert details.get("status") == 503
    assert details.get("retry_after_source") == "header"
    assert details.get("retry_after") == pytest.approx(30.0, abs=1e-3)

    cache = pubchem._ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
    entry = cache.get(pubchem._build_cache_key("GET", url))
    assert entry is not None
    assert entry.outcome == "server_error"
    assert entry.details == {
        "reason": "server_error",
        "status": 503,
        "retry_after": pytest.approx(30.0, abs=1e-3),
        "retry_after_source": "header",
        "cache": True,
    }

    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_UNTIL", None)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_DETAILS", None)

    second_result = pubchem.make_request(url, cfg)

    assert second_result is None
    assert len(session.calls) == 1
    assert limiter.acquires == 1
    assert sleep_calls == []
    outcome, details = pubchem.last_request_outcome()
    assert outcome == "server_error"
    assert details == {
        "reason": "server_error",
        "retry_after": pytest.approx(30.0, abs=1e-3),
        "retry_after_source": "header",
        "http_status": 503,
        "cache": True,
    }


def test_make_request__short_circuits_during_retry_after(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()
    cfg.retries = 0

    base = cfg.base.rstrip("/")
    url_first = (
        f"{base}/compound/cid/149790/property/"
        "MolecularFormula,IUPACName,SMILES,ConnectivitySMILES,InChI,InChIKey/JSON"
    )
    url_second = (
        f"{base}/compound/cid/20353/property/"
        "MolecularFormula,IUPACName,SMILES,ConnectivitySMILES,InChI,InChIKey/JSON"
    )

    time_ctrl = _TimeController()
    monkeypatch.setattr(pubchem, "monotonic", time_ctrl.monotonic)
    monkeypatch.setattr(pubchem, "sleep", time_ctrl.sleep)

    limiter = _DummyLimiter()
    monkeypatch.setattr(pubchem, "get_limiter", lambda *_args, **_kwargs: limiter)

    success_payload = {
        "PropertyTable": {
            "Properties": [
                {
                    "MolecularFormula": "C9H8O4",
                    "IUPACName": "2-acetyloxybenzoic acid",
                }
            ]
        }
    }

    responses = [
        _DummyResponse(503, {"Retry-After": "30"}),
        _DummyResponse(200, {}, success_payload),
    ]
    index = {"value": 0}

    def _response() -> _DummyResponse:
        value = responses[index["value"]]
        index["value"] += 1
        return value

    session = _DummySession(_response)
    monkeypatch.setattr(pubchem, "get_session", lambda *_args, **_kwargs: session)
    monkeypatch.setattr(pubchem, "_CACHE", None)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_UNTIL", None)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_DETAILS", None)

    first = pubchem.make_request(url_first, cfg)

    assert first is None
    assert session.calls == [
        (
            "GET",
            url_first,
            {"timeout": (cfg.timeout_connect, cfg.timeout_read)},
        )
    ]
    assert limiter.acquires == 1
    outcome, details = pubchem.last_request_outcome()
    assert outcome == "server_error"
    assert details is not None
    assert details.get("retry_after_source") == "header"

    second = pubchem.make_request(url_second, cfg)

    assert second is None
    assert session.calls == [
        (
            "GET",
            url_first,
            {"timeout": (cfg.timeout_connect, cfg.timeout_read)},
        )
    ]
    assert limiter.acquires == 1
    outcome, details = pubchem.last_request_outcome()
    assert outcome == "server_error"
    assert details is not None
    assert details.get("retry_after_source") == "header"
    assert details.get("retry_after") == pytest.approx(30.0)
    assert details.get("cache") is True

    time_ctrl.advance(31.0)

    third = pubchem.make_request(url_second, cfg)

    assert third == success_payload
    assert session.calls == [
        (
            "GET",
            url_first,
            {"timeout": (cfg.timeout_connect, cfg.timeout_read)},
        ),
        (
            "GET",
            url_second,
            {"timeout": (cfg.timeout_connect, cfg.timeout_read)},
        ),
    ]
    assert limiter.acquires == 2
    outcome, details = pubchem.last_request_outcome()
    assert outcome == "hit"
    assert details == {"status": 200}

def test_make_request__invalid_identifier_cached(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
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
    outcome, details = pubchem.last_request_outcome()
    assert outcome == "invalid_identifier"
    assert details == {"reason": "invalid_identifier", "status": 400}

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
    outcome, details = pubchem.last_request_outcome()
    assert outcome == "invalid_identifier"
    assert details == {
        "cache": True,
        "reason": "invalid_identifier",
        "http_status": 400,
    }


def test_make_request__retry_after_honours_grace(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()
    cfg.retries = 3
    cfg.timeout_seconds = 10

    url = (
        f"{cfg.base.rstrip('/')}/compound/cid/64972/property/"
        "MolecularFormula,IUPACName,IsomericSMILES,CanonicalSMILES,InChI,InChIKey/JSON"
    )

    sleep_calls: list[float] = []
    monkeypatch.setattr(pubchem, "sleep", lambda seconds: sleep_calls.append(seconds))
    limiter = _DummyLimiter()
    monkeypatch.setattr(pubchem, "get_limiter", lambda *_args, **_kwargs: limiter)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_UNTIL", None)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_DETAILS", None)

    def _response() -> _DummyResponse:
        return _DummyResponse(503, {"Retry-After": "30"})

    session = _DummySession(_response)
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
    assert sleep_calls == []
    assert limiter.acquires == 1
    outcome, details = pubchem.last_request_outcome()
    assert outcome == "server_error"
    assert details == {
        "reason": "server_error",
        "status": 503,
        "retry_after": pytest.approx(30.0, abs=1e-3),
        "retry_after_source": "header",
        "cache": True,
    }

    cache = pubchem._ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
    entry = cache.get(pubchem._build_cache_key("GET", url))
    assert entry is not None
    assert entry.outcome == "server_error"
    assert entry.details == {
        "reason": "server_error",
        "status": 503,
        "retry_after": pytest.approx(30.0, abs=1e-3),
        "retry_after_source": "header",
        "cache": True,
    }


def test_make_request__retry_after_grace_disabled_causes_timeout(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()
    cfg.retries = 3
    cfg.timeout_seconds = 10
    cfg.retry_after_grace_seconds = 0

    url = (
        f"{cfg.base.rstrip('/')}/compound/cid/64972/property/"
        "MolecularFormula,IUPACName,IsomericSMILES,CanonicalSMILES,InChI,InChIKey/JSON"
    )

    sleep_calls: list[float] = []
    monkeypatch.setattr(pubchem, "sleep", lambda seconds: sleep_calls.append(seconds))
    limiter = _DummyLimiter()
    monkeypatch.setattr(pubchem, "get_limiter", lambda *_args, **_kwargs: limiter)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_UNTIL", None)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_DETAILS", None)

    def _response() -> _DummyResponse:
        return _DummyResponse(503, {"Retry-After": "30"})

    session = _DummySession(_response)
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
    assert sleep_calls == []
    assert limiter.acquires == 1
    outcome, details = pubchem.last_request_outcome()
    assert outcome == "timeout"
    assert details is not None
    assert details.get("timeout_reason") == "retry_after_exceeds_deadline"

    cache = pubchem._ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
    entry = cache.get(pubchem._build_cache_key("GET", url))
    assert entry is not None
    assert entry.outcome == "timeout"
    assert entry.details is not None
    details = dict(entry.details)
    timeout_retry_after = details.pop("timeout_retry_after", None)
    timeout_flag = details.pop("timeout", None)
    timeout_stored_at = details.pop("timeout_stored_at", None)
    assert timeout_flag is True
    assert isinstance(timeout_stored_at, float)
    assert timeout_retry_after == pytest.approx(30.0)
    assert details == {
        "reason": "server_error",
        "status": 503,
        "retry_after": 30.0,
        "retry_after_source": "header",
        "timeout_reason": "retry_after_exceeds_deadline",
    }


def test_make_request__applies_jitter_without_retry_after(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()
    cfg.retries = 1
    cfg.retry_jitter_seconds = 0.35
    cfg.retry_jitter_seed = 17

    url = (
        f"{cfg.base.rstrip('/')}/compound/cid/64972/property/"
        "MolecularFormula,IUPACName,IsomericSMILES,CanonicalSMILES,InChI,InChIKey/JSON"
    )

    limiter = _DummyLimiter()
    monkeypatch.setattr(pubchem, "get_limiter", lambda *_args, **_kwargs: limiter)
    monkeypatch.setattr(pubchem, "_CACHE", None)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_UNTIL", None)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_DETAILS", None)
    monkeypatch.setattr(pubchem, "_JITTER_GENERATORS", {})

    sleep_calls: list[float] = []
    monkeypatch.setattr(pubchem, "sleep", lambda seconds: sleep_calls.append(float(seconds)))

    session = _DummySession(lambda: _DummyResponse(503, {}))
    monkeypatch.setattr(pubchem, "get_session", lambda *_args, **_kwargs: session)

    result = pubchem.make_request(url, cfg)

    assert result is None
    assert limiter.acquires == 1
    assert len(sleep_calls) == 1
    expected_jitter = random.Random(cfg.retry_jitter_seed).uniform(
        0.0, cfg.retry_jitter_seconds
    )
    expected_delay = cfg.backoff_initial_seconds + expected_jitter
    assert sleep_calls[0] == pytest.approx(expected_delay, rel=1e-9)


def test_make_request__timeout_cache_uses_config_backoff(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()
    cfg.retries = 3
    cfg.timeout_seconds = 5
    cfg.backoff_initial_seconds = 10.0

    url = (
        f"{cfg.base.rstrip('/')}/compound/cid/64972/property/"
        "MolecularFormula,IUPACName,IsomericSMILES,CanonicalSMILES,InChI,InChIKey/JSON"
    )

    sleep_calls: list[float] = []
    monkeypatch.setattr(
        pubchem, "sleep", lambda seconds: sleep_calls.append(float(seconds))
    )
    limiter = _DummyLimiter()
    monkeypatch.setattr(pubchem, "get_limiter", lambda *_args, **_kwargs: limiter)

    session = _DummySession(lambda: _DummyResponse(503, {}))
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
    assert sleep_calls == []
    assert limiter.acquires == 1
    cache = pubchem._ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
    entry = cache.get(pubchem._build_cache_key("GET", url))
    assert entry is not None
    assert entry.outcome == "timeout"
    assert entry.details is not None
    details = dict(entry.details)
    timeout_retry_after = details.pop("timeout_retry_after", None)
    timeout_flag = details.pop("timeout", None)
    timeout_stored_at = details.pop("timeout_stored_at", None)
    assert timeout_flag is True
    assert isinstance(timeout_stored_at, float)
    assert timeout_retry_after == pytest.approx(cfg.timeout_seconds)
    assert details == {
        "reason": "server_error",
        "status": 503,
        "retry_after": pytest.approx(cfg.backoff_initial_seconds),
        "retry_after_source": "backoff",
        "timeout_reason": "retry_after_exceeds_deadline",
    }


    cache = pubchem._ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
    entry = cache.get(pubchem._build_cache_key("GET", url))
    assert entry is not None
    assert entry.outcome == "timeout"
    assert entry.details is not None
    details = dict(entry.details)
    timeout_retry_after = details.pop("timeout_retry_after", None)
    timeout_flag = details.pop("timeout", None)
    timeout_stored_at = details.pop("timeout_stored_at", None)
    assert timeout_flag is True
    assert isinstance(timeout_stored_at, float)
    assert timeout_retry_after == pytest.approx(cfg.timeout_seconds)
    assert details == {
        "reason": "server_error",
        "status": 503,
        "retry_after": pytest.approx(cfg.backoff_initial_seconds),
        "retry_after_source": "backoff",
        "timeout_reason": "retry_after_exceeds_deadline",
    }


def test_get_cid_from_smiles__uses_post_for_stereochemistry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()
    smiles = r"C/C=C\C"

    limiter = _DummyLimiter()
    monkeypatch.setattr(pubchem, "get_limiter", lambda *_args, **_kwargs: limiter)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_UNTIL", None)
    monkeypatch.setattr(pubchem, "_SERVICE_UNAVAILABLE_DETAILS", None)

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


def test_get_cid__falls_back_to_pug_when_rdf_fails(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()
    calls: list[str] = []

    def _fake_make_request(url: str, _cfg: PubChemCfg, **_kwargs: object) -> dict[str, Any] | None:
        calls.append(url)
        if "/rdf/query" in url:
            return None
        assert url.endswith(
            f"/compound/name/{pubchem.url_encode('Example')}/cids/JSON"
        )
        return {"IdentifierList": {"CID": ["10", "2", "10"]}}

    monkeypatch.setattr(pubchem, "make_request", _fake_make_request)

    result = pubchem.get_cid("Example", cfg)

    assert result == "10|2"
    assert any("/rdf/query" in call for call in calls)
    assert any("/compound/name/" in call for call in calls)


def test_get_all_cid__falls_back_to_pug_when_rdf_fails(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()
    calls: list[str] = []

    def _fake_make_request(url: str, _cfg: PubChemCfg, **_kwargs: object) -> dict[str, Any] | None:
        calls.append(url)
        if "/rdf/query" in url:
            return None
        assert url.endswith(
            f"/compound/name/{pubchem.url_encode('Sample')}/cids/JSON?name_type=word"
        )
        return {"IdentifierList": {"CID": ["77", "42"]}}

    monkeypatch.setattr(pubchem, "make_request", _fake_make_request)

    result = pubchem.get_all_cid("Sample", cfg)

    assert result == "42|77"
    assert any("/rdf/query" in call for call in calls)
    assert any("name_type=word" in call for call in calls)


def test_get_cid_from_smiles__raises_on_service_unavailable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()
    outcome_details = {"status": 503, "reason": "server_error"}

    monkeypatch.setattr(pubchem, "make_request", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        pubchem, "last_request_outcome", lambda: ("server_error", outcome_details)
    )

    with pytest.raises(pubchem.PubChemServiceUnavailable) as excinfo:
        pubchem.get_cid_from_smiles("CCO", cfg)

    assert excinfo.value.outcome == "server_error"
    assert excinfo.value.status == 503
    assert excinfo.value.details == outcome_details


def test_get_cid_from_smiles__returns_none_on_not_found(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = PubChemCfg()

    monkeypatch.setattr(pubchem, "make_request", lambda *args, **kwargs: None)
    monkeypatch.setattr(pubchem, "last_request_outcome", lambda: ("not_found", {}))

    assert pubchem.get_cid_from_smiles("CCO", cfg) is None
