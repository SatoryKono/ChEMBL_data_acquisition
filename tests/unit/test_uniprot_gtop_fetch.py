from __future__ import annotations

import json
from collections.abc import Callable, Iterator
from types import SimpleNamespace

import pytest
import requests

from library.config import IupharCfg
from library.integration import uniprot_library as uniprot


@pytest.fixture()
def reset_gtop_caches() -> Iterator[None]:
    """Ensure Guide-to-Pharmacology fetch caches are cleared for each test."""

    uniprot._GTOP_JSON_FAILURE_CACHE.clear()
    uniprot._GTOP_NON_JSON_CONTENT_TYPE_CACHE.clear()
    uniprot._GTOP_SKIPPED_FAILURE_LOG.clear()
    yield
    uniprot._GTOP_JSON_FAILURE_CACHE.clear()
    uniprot._GTOP_NON_JSON_CONTENT_TYPE_CACHE.clear()
    uniprot._GTOP_SKIPPED_FAILURE_LOG.clear()


class _DummyResponse:
    def __init__(
        self,
        payload: object,
        content_type: str,
        *,
        json_exc: Exception | None = None,
    ):
        self._payload = payload
        self._json_exc = json_exc
        self.headers = {"Content-Type": content_type}
        self.status_code = 200

    def raise_for_status(self) -> None:  # pragma: no cover - always OK in tests
        return None

    @property
    def content(self) -> bytes:
        if isinstance(self._payload, (bytes, bytearray)):
            return bytes(self._payload)
        return json.dumps(self._payload).encode("utf-8")

    def json(self) -> object:
        if self._json_exc is not None:
            raise self._json_exc
        return self._payload

    def __enter__(self) -> "_DummyResponse":  # pragma: no cover - context manager boilerplate
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:  # pragma: no cover - context manager boilerplate
        return False


class _DummySession:
    def __init__(self, response_factory: Callable[[], _DummyResponse]):
        self._factory = response_factory
        self.calls: int = 0
        self.last_url: str | None = None
        self.last_timeout: tuple[float, float] | None = None

    def get(self, url: str, timeout: tuple[float, float]) -> _DummyResponse:
        self.calls += 1
        self.last_url = url
        self.last_timeout = timeout
        return self._factory()


class _FailingSession:
    def __init__(self, exc: Exception):
        self._exc = exc
        self.calls = 0
        self.last_url: str | None = None
        self.last_timeout: tuple[float, float] | None = None

    def get(self, url: str, timeout: tuple[float, float]) -> _DummyResponse:
        self.calls += 1
        self.last_url = url
        self.last_timeout = timeout
        raise self._exc


def _patch_dependencies(monkeypatch: pytest.MonkeyPatch, session: _DummySession) -> None:
    monkeypatch.setattr(uniprot, "get_uniprot_session", lambda: session)
    monkeypatch.setattr(
        uniprot,
        "get_limiter",
        lambda name, rps, burst=None: SimpleNamespace(acquire=lambda: None),
    )


@pytest.mark.unit
def test_fetch_gtop_endpoint__accepts_text_plain_json(
    monkeypatch: pytest.MonkeyPatch,
    reset_gtop_caches: None,
) -> None:
    session = _DummySession(
        lambda: _DummyResponse({"value": 1}, "text/plain; charset=utf-8")
    )
    _patch_dependencies(monkeypatch, session)

    events: list[tuple[str, dict[str, object]]] = []

    def capture(event: str, *args: object, **kwargs: object) -> None:
        events.append((event, dict(kwargs)))

    monkeypatch.setattr(uniprot.logger, "warning", capture)

    cfg = IupharCfg()

    result = uniprot._fetch_gtop_endpoint("GTP1", "function", cfg=cfg)
    assert result == {"value": 1}
    assert events == [
        (
            "gtop_non_json_response",
            {
                "gtop_id": "GTP1",
                "endpoint": "function",
                "content_type": "text/plain; charset=utf-8",
            },
        )
    ]
    assert ("GTP1", "function") not in uniprot._GTOP_JSON_FAILURE_CACHE
    assert ("GTP1", "function") in uniprot._GTOP_NON_JSON_CONTENT_TYPE_CACHE
    assert session.calls == 1

    result_again = uniprot._fetch_gtop_endpoint("GTP1", "function", cfg=cfg)
    assert result_again == {"value": 1}
    assert events == [
        (
            "gtop_non_json_response",
            {
                "gtop_id": "GTP1",
                "endpoint": "function",
                "content_type": "text/plain; charset=utf-8",
            },
        )
    ]
    assert session.calls == 2


@pytest.mark.unit
def test_fetch_gtop_endpoint__caches_request_exception(
    monkeypatch: pytest.MonkeyPatch,
    reset_gtop_caches: None,
) -> None:
    session = _FailingSession(requests.RequestException("boom"))
    _patch_dependencies(monkeypatch, session)

    warning_events: list[tuple[str, dict[str, object]]] = []

    def capture_warning(event: str, *args: object, **kwargs: object) -> None:
        warning_events.append((event, dict(kwargs)))

    monkeypatch.setattr(uniprot.logger, "warning", capture_warning)

    cfg = IupharCfg()

    result = uniprot._fetch_gtop_endpoint("GTP4", "function", cfg=cfg)
    assert result is None
    assert warning_events == [
        (
            "gtop_request_failed",
            {"gtop_id": "GTP4", "endpoint": "function", "error": "boom"},
        )
    ]
    assert session.calls == 1
    assert ("GTP4", "function") in uniprot._GTOP_JSON_FAILURE_CACHE

    warning_events.clear()

    second = uniprot._fetch_gtop_endpoint("GTP4", "function", cfg=cfg)
    assert second is None
    assert warning_events == []
    assert session.calls == 1


@pytest.mark.unit
def test_fetch_gtop_endpoint__handles_empty_body(
    monkeypatch: pytest.MonkeyPatch,
    reset_gtop_caches: None,
) -> None:
    response = _DummyResponse(b"", "application/json")
    response.headers["Content-Length"] = "0"
    session = _DummySession(lambda: response)
    _patch_dependencies(monkeypatch, session)

    info_events: list[tuple[str, dict[str, object]]] = []
    warning_events: list[tuple[str, dict[str, object]]] = []

    def capture_info(event: str, *args: object, **kwargs: object) -> None:
        info_events.append((event, dict(kwargs)))

    def capture_warning(event: str, *args: object, **kwargs: object) -> None:
        warning_events.append((event, dict(kwargs)))

    monkeypatch.setattr(uniprot.logger, "info", capture_info)
    monkeypatch.setattr(uniprot.logger, "warning", capture_warning)

    cfg = IupharCfg()

    result = uniprot._fetch_gtop_endpoint("GTP3", "naturalLigands", cfg=cfg)
    assert result == []
    assert info_events == [
        (
            "gtop_empty_response",
            {
                "gtop_id": "GTP3",
                "endpoint": "naturalLigands",
                "content_type": "application/json",
            },
        )
    ]
    assert warning_events == []
    assert ("GTP3", "naturalLigands") not in uniprot._GTOP_JSON_FAILURE_CACHE
    assert session.calls == 1

    second = uniprot._fetch_gtop_endpoint("GTP3", "naturalLigands", cfg=cfg)
    assert second == []
    assert info_events == [
        (
            "gtop_empty_response",
            {
                "gtop_id": "GTP3",
                "endpoint": "naturalLigands",
                "content_type": "application/json",
            },
        ),
        (
            "gtop_empty_response",
            {
                "gtop_id": "GTP3",
                "endpoint": "naturalLigands",
                "content_type": "application/json",
            },
        ),
    ]
    assert warning_events == []
    assert ("GTP3", "naturalLigands") not in uniprot._GTOP_JSON_FAILURE_CACHE
    assert session.calls == 2


@pytest.mark.unit
def test_fetch_gtop_endpoint__records_decode_failure(
    monkeypatch: pytest.MonkeyPatch,
    reset_gtop_caches: None,
) -> None:
    json_error = json.JSONDecodeError("invalid", "{", 0)
    session = _DummySession(
        lambda: _DummyResponse({"invalid": True}, "application/json", json_exc=json_error)
    )
    _patch_dependencies(monkeypatch, session)

    events: list[tuple[str, dict[str, object]]] = []

    def capture(event: str, *args: object, **kwargs: object) -> None:
        events.append((event, dict(kwargs)))

    monkeypatch.setattr(uniprot.logger, "warning", capture)

    cfg = IupharCfg()

    result = uniprot._fetch_gtop_endpoint("GTP2", "interactions", cfg=cfg)
    assert result is None
    assert events == [
        (
            "gtop_json_decode_failed",
            {
                "gtop_id": "GTP2",
                "endpoint": "interactions",
                "error": "invalid: line 1 column 1 (char 0)",
                "content_type": "application/json",
            },
        )
    ]
    assert ("GTP2", "interactions") in uniprot._GTOP_JSON_FAILURE_CACHE
    assert session.calls == 1

    second = uniprot._fetch_gtop_endpoint("GTP2", "interactions", cfg=cfg)
    assert second is None
    assert events == [
        (
            "gtop_json_decode_failed",
            {
                "gtop_id": "GTP2",
                "endpoint": "interactions",
                "error": "invalid: line 1 column 1 (char 0)",
                "content_type": "application/json",
            },
        )
    ]
    assert session.calls == 1
