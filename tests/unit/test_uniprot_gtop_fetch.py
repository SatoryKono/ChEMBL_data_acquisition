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
    uniprot._reset_gtop_circuit_state()
    yield
    uniprot._GTOP_JSON_FAILURE_CACHE.clear()
    uniprot._GTOP_NON_JSON_CONTENT_TYPE_CACHE.clear()
    uniprot._GTOP_SKIPPED_FAILURE_LOG.clear()
    uniprot._reset_gtop_circuit_state()


class _DummyResponse:
    def __init__(
        self,
        payload: object,
        content_type: str,
        *,
        json_exc: Exception | None = None,
        status_code: int = 200,
        text: str | None = None,
    ):
        self._payload = payload
        self._json_exc = json_exc
        self.headers = {"Content-Type": content_type}
        self.status_code = status_code
        self._text_override = text

    def raise_for_status(self) -> None:  # pragma: no cover - always OK in tests
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code} error", response=self)

    @property
    def content(self) -> bytes:
        if isinstance(self._payload, bytes | bytearray):
            return bytes(self._payload)
        if self._text_override is not None:
            return self._text_override.encode("utf-8")
        return json.dumps(self._payload).encode("utf-8")

    @property
    def text(self) -> str:
        if self._text_override is not None:
            return self._text_override
        if isinstance(self._payload, str):
            return self._payload
        return json.dumps(self._payload)

    def json(self) -> object:
        if self._json_exc is not None:
            raise self._json_exc
        return self._payload

    def __enter__(
        self,
    ) -> _DummyResponse:  # pragma: no cover - context manager boilerplate
        return self

    def __exit__(
        self, exc_type, exc, tb
    ) -> bool:  # pragma: no cover - context manager boilerplate
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


def _patch_dependencies(
    monkeypatch: pytest.MonkeyPatch, session: _DummySession
) -> None:
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
            {
                "gtop_id": "GTP4",
                "endpoint": "function",
                "error": "boom",
                "status_code": None,
            },
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
def test_fetch_gtop_endpoint__skips_when_disabled(
    monkeypatch: pytest.MonkeyPatch,
    reset_gtop_caches: None,
) -> None:
    session = _DummySession(lambda: _DummyResponse({}, "application/json"))
    _patch_dependencies(monkeypatch, session)

    debug_events: list[tuple[str, dict[str, object]]] = []

    def capture_debug(event: str, *args: object, **kwargs: object) -> None:
        debug_events.append((event, dict(kwargs)))

    monkeypatch.setattr(uniprot.logger, "debug", capture_debug)

    cfg = IupharCfg(enable=False)

    result = uniprot._fetch_gtop_endpoint("GTP9", "function", cfg=cfg)

    assert result == []
    assert session.calls == 0
    assert (
        "gtop_fetch_disabled",
        {"gtop_id": "GTP9", "endpoint": "function"},
    ) in debug_events


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
        lambda: _DummyResponse(
            {"invalid": True}, "application/json", json_exc=json_error
        )
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


@pytest.mark.unit
def test_update_gtop_metadata__disabled_skips_fetch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    result = {
        "GuidetoPHARMACOLOGY": "GTP9",
        "gtop_natural_ligands_n": "",
        "gtop_interactions_n": "",
        "gtop_function_text_short": "",
    }

    def fail(*args: object, **kwargs: object) -> None:  # pragma: no cover - fails test
        raise AssertionError("_fetch_gtop_endpoint should not be called when disabled")

    monkeypatch.setattr(uniprot, "_fetch_gtop_endpoint", fail)

    cfg = IupharCfg(enable=False)

    uniprot._update_gtop_metadata(result, cfg=cfg)

    assert result["gtop_natural_ligands_n"] == ""
    assert result["gtop_interactions_n"] == ""
    assert result["gtop_function_text_short"] == ""


def test_fetch_gtop_endpoint__handles_404_missing_target(
    monkeypatch: pytest.MonkeyPatch,
    reset_gtop_caches: None,
) -> None:
    session = _DummySession(
        lambda: _DummyResponse({}, "application/json", status_code=404, text="{}")
    )
    _patch_dependencies(monkeypatch, session)

    info_events: list[tuple[str, dict[str, object]]] = []
    warning_events: list[tuple[str, dict[str, object]]] = []

    monkeypatch.setattr(
        uniprot.logger, "info", lambda event, **kw: info_events.append((event, kw))
    )
    monkeypatch.setattr(
        uniprot.logger,
        "warning",
        lambda event, **kw: warning_events.append((event, kw)),
    )

    cfg = IupharCfg()

    result = uniprot._fetch_gtop_endpoint("GTP5", "interactions", cfg=cfg)
    assert result is None
    assert info_events == [
        (
            "gtop_endpoint_missing",
            {"gtop_id": "GTP5", "endpoint": "interactions", "status_code": 404},
        )
    ]
    assert warning_events == []
    assert ("GTP5", "interactions") in uniprot._GTOP_JSON_FAILURE_CACHE
    assert session.calls == 1

    second = uniprot._fetch_gtop_endpoint("GTP5", "interactions", cfg=cfg)
    assert second is None
    assert session.calls == 1


class _FakeMonotonic:
    def __init__(self) -> None:
        self.value = 0.0

    def __call__(self) -> float:
        return self.value

    def advance(self, delta: float) -> None:
        self.value += delta


@pytest.mark.unit
def test_fetch_gtop_endpoint__circuit_breaker_limits_retries(
    monkeypatch: pytest.MonkeyPatch,
    reset_gtop_caches: None,
) -> None:
    response = requests.Response()
    response.status_code = 503
    response.url = "https://example"
    exc = requests.HTTPError("503 Server Error", response=response)
    session = _FailingSession(exc)
    _patch_dependencies(monkeypatch, session)

    clock = _FakeMonotonic()
    monkeypatch.setattr(uniprot.time, "monotonic", clock)
    monkeypatch.setattr(uniprot, "_GTOP_CIRCUIT_FAILURE_THRESHOLD", 2)
    monkeypatch.setattr(uniprot, "_GTOP_CIRCUIT_HOLDOFF_SECONDS", 10.0)

    warning_events: list[tuple[str, dict[str, object]]] = []
    info_events: list[tuple[str, dict[str, object]]] = []

    def capture_warning(event: str, **kwargs: object) -> None:
        warning_events.append((event, dict(kwargs)))

    def capture_info(event: str, **kwargs: object) -> None:
        info_events.append((event, dict(kwargs)))

    monkeypatch.setattr(uniprot.logger, "warning", capture_warning)
    monkeypatch.setattr(uniprot.logger, "info", capture_info)

    cfg = IupharCfg()

    first_id = "GTP6A"
    assert uniprot._fetch_gtop_endpoint(first_id, "function", cfg=cfg) is None
    assert session.calls == 1
    assert warning_events == [
        (
            "gtop_request_failed",
            {
                "gtop_id": first_id,
                "endpoint": "function",
                "error": "503 Server Error",
                "status_code": 503,
            },
        )
    ]

    warning_events.clear()

    second_id = "GTP6B"
    assert uniprot._fetch_gtop_endpoint(second_id, "function", cfg=cfg) is None
    assert session.calls == 2
    assert warning_events == [
        (
            "gtop_request_failed",
            {
                "gtop_id": second_id,
                "endpoint": "function",
                "error": "503 Server Error",
                "status_code": 503,
            },
        ),
        (
            "gtop_circuit_opened",
            {
                "gtop_id": second_id,
                "endpoint": "function",
                "retry_after": 10.0,
                "failure_threshold": 2,
            },
        ),
    ]

    clock.advance(1.0)
    warning_events.clear()

    third_id = "GTP6C"
    assert uniprot._fetch_gtop_endpoint(third_id, "function", cfg=cfg) is None
    assert session.calls == 2
    assert warning_events == [
        (
            "gtop_circuit_open_skip",
            {"gtop_id": third_id, "endpoint": "function", "retry_after": 9.0},
        )
    ]

    warning_events.clear()
    info_events.clear()

    clock.advance(10.0)

    fourth_id = "GTP6D"
    assert uniprot._fetch_gtop_endpoint(fourth_id, "function", cfg=cfg) is None
    assert session.calls == 3
    assert info_events == [("gtop_circuit_reset", {"downtime": 10.0})]
    assert warning_events == [
        (
            "gtop_request_failed",
            {
                "gtop_id": fourth_id,
                "endpoint": "function",
                "error": "503 Server Error",
                "status_code": 503,
            },
        )
    ]

    warning_events.clear()

    fifth_id = "GTP6E"
    assert uniprot._fetch_gtop_endpoint(fifth_id, "function", cfg=cfg) is None
    assert session.calls == 4
    assert warning_events == [
        (
            "gtop_request_failed",
            {
                "gtop_id": fifth_id,
                "endpoint": "function",
                "error": "503 Server Error",
                "status_code": 503,
            },
        ),
        (
            "gtop_circuit_opened",
            {
                "gtop_id": fifth_id,
                "endpoint": "function",
                "retry_after": 10.0,
                "failure_threshold": 2,
            },
        ),
    ]

    clock.advance(1.0)
    warning_events.clear()

    sixth_id = "GTP6F"
    assert uniprot._fetch_gtop_endpoint(sixth_id, "function", cfg=cfg) is None
    assert session.calls == 4
    assert warning_events == [
        (
            "gtop_circuit_open_skip",
            {"gtop_id": sixth_id, "endpoint": "function", "retry_after": 9.0},
        )
    ]


@pytest.mark.unit
def test_fetch_gtop_endpoint__circuit_breaker_counts_per_endpoint(
    monkeypatch: pytest.MonkeyPatch,
    reset_gtop_caches: None,
) -> None:
    session: _DummySession

    def response_factory() -> _DummyResponse:
        assert session.last_url is not None
        if session.last_url.endswith("/naturalLigands"):
            return _DummyResponse([], "application/json")
        if session.last_url.endswith("/interactions"):
            return _DummyResponse([], "application/json")
        if session.last_url.endswith("/function"):
            return _DummyResponse({}, "application/json", status_code=500)
        raise AssertionError(f"unexpected URL {session.last_url}")

    session = _DummySession(response_factory)
    _patch_dependencies(monkeypatch, session)

    clock = _FakeMonotonic()
    monkeypatch.setattr(uniprot.time, "monotonic", clock)
    monkeypatch.setattr(uniprot, "_GTOP_CIRCUIT_FAILURE_THRESHOLD", 2)
    monkeypatch.setattr(uniprot, "_GTOP_CIRCUIT_HOLDOFF_SECONDS", 10.0)

    warning_events: list[tuple[str, dict[str, object]]] = []
    monkeypatch.setattr(
        uniprot.logger,
        "warning",
        lambda event, **kwargs: warning_events.append((event, dict(kwargs))),
    )

    cfg = IupharCfg()

    # First ID: successes for other endpoints should not affect the function breaker
    assert uniprot._fetch_gtop_endpoint("GTP7A", "naturalLigands", cfg=cfg) == []
    assert uniprot._fetch_gtop_endpoint("GTP7A", "interactions", cfg=cfg) == []
    assert uniprot._fetch_gtop_endpoint("GTP7A", "function", cfg=cfg) is None
    assert warning_events == [
        (
            "gtop_request_failed",
            {
                "gtop_id": "GTP7A",
                "endpoint": "function",
                "error": "500 error",
                "status_code": 500,
            },
        )
    ]

    warning_events.clear()

    # Second ID: after another failure for the same endpoint the circuit opens
    assert uniprot._fetch_gtop_endpoint("GTP7B", "naturalLigands", cfg=cfg) == []
    assert uniprot._fetch_gtop_endpoint("GTP7B", "interactions", cfg=cfg) == []
    assert uniprot._fetch_gtop_endpoint("GTP7B", "function", cfg=cfg) is None
    assert warning_events == [
        (
            "gtop_request_failed",
            {
                "gtop_id": "GTP7B",
                "endpoint": "function",
                "error": "500 error",
                "status_code": 500,
            },
        ),
        (
            "gtop_circuit_opened",
            {
                "gtop_id": "GTP7B",
                "endpoint": "function",
                "retry_after": 10.0,
                "failure_threshold": 2,
            },
        ),
    ]

    warning_events.clear()
    clock.advance(1.0)

    # Third ID: circuit remains open for the function endpoint and skips the call
    assert uniprot._fetch_gtop_endpoint("GTP7C", "naturalLigands", cfg=cfg) == []
    assert uniprot._fetch_gtop_endpoint("GTP7C", "interactions", cfg=cfg) == []
    assert uniprot._fetch_gtop_endpoint("GTP7C", "function", cfg=cfg) is None
    assert warning_events == [
        (
            "gtop_circuit_open_skip",
            {"gtop_id": "GTP7C", "endpoint": "function", "retry_after": 9.0},
        )
    ]

    assert session.calls == 8
