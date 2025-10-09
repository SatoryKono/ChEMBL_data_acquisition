from __future__ import annotations

import pytest
import requests
from urllib3.util import retry as urllib3_retry

from library.clients import pubmed
from library.config.models import ApiCfg, PubMedCfg, RetryCfg
from library.config.runtime import session_with_retry


class _DummyResponse:
    def __init__(self, text: str = "<PubmedArticleSet></PubmedArticleSet>") -> None:
        self.status_code = 200
        self.text = text
        self.headers: dict[str, str] = {}

    def json(self) -> dict[str, str]:  # pragma: no cover - JSON not expected here
        raise AssertionError("JSON payloads are not used in PubMed XML fetches")

    def __enter__(self) -> "_DummyResponse":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: object | None,
    ) -> bool:
        return False


class _DummySession:
    def __init__(self) -> None:
        self.calls: list[tuple[str, float | tuple[float, float], dict[str, object]]] = []

    def get(self, url: str, *, timeout: float | tuple[float, float], **kwargs: object) -> _DummyResponse:
        self.calls.append((url, timeout, kwargs))
        return _DummyResponse()


class _SequencedResponse:
    def __init__(
        self,
        status_code: int,
        text: str = "",
        headers: dict[str, str] | None = None,
    ) -> None:
        self.status_code = status_code
        self.text = text
        self.headers = headers or {}

    def __enter__(self) -> "_SequencedResponse":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: object | None,
    ) -> bool:
        return False


class _SequencedSession:
    def __init__(self, responses: list[_SequencedResponse]) -> None:
        self._responses = responses
        self.calls: list[tuple[str, float | tuple[float, float], dict[str, object]]] = []

    def get(
        self,
        url: str,
        *,
        timeout: float | tuple[float, float],
        **kwargs: object,
    ) -> _SequencedResponse:
        if not self._responses:
            raise AssertionError("No responses left in sequence")
        self.calls.append((url, timeout, kwargs))
        return self._responses.pop(0)


@pytest.mark.unit
def test_retry_delay__respects_backoff_cap() -> None:
    retry_cfg = RetryCfg(max_attempts=4, backoff_factor=2.0, backoff_cap=3.0)

    delay = pubmed._retry_delay(
        2,
        0.25,
        retry_cfg,
        timeout=None,
        jitter=lambda _: 0.5,
    )

    assert delay == pytest.approx(3.0)


@pytest.mark.unit
def test_retry_delay__deterministic_jitter() -> None:
    retry_cfg = RetryCfg(max_attempts=4, backoff_factor=1.0, backoff_cap=None, jitter_seed=7)
    jitter_one = retry_cfg.build_jitter()
    jitter_two = RetryCfg(max_attempts=4, backoff_factor=1.0, backoff_cap=None, jitter_seed=7).build_jitter()

    assert jitter_one is not None
    assert jitter_two is not None

    base_delay = 0.25
    delays_one = [
        pubmed._retry_delay(attempt, base_delay, retry_cfg, timeout=None, jitter=jitter_one)
        for attempt in range(1, 4)
    ]
    delays_two = [
        pubmed._retry_delay(attempt, base_delay, retry_cfg, timeout=None, jitter=jitter_two)
        for attempt in range(1, 4)
    ]

    assert delays_one == delays_two


@pytest.mark.unit
def test_session_with_retry__uses_backoff_cap() -> None:
    api_cfg = ApiCfg()
    retry_cfg = RetryCfg(max_attempts=4, backoff_factor=1.0, backoff_cap=7.5)

    with session_with_retry(api_cfg, retry_cfg) as session:
        adapter = session.get_adapter("https://")
        max_retries = adapter.max_retries

    assert isinstance(max_retries, urllib3_retry.Retry)
    assert max_retries.backoff_max == pytest.approx(7.5)


@pytest.mark.unit
def test_fetch_pubmed_batch__includes_contact_parameters() -> None:
    session = _DummySession()
    cfg = PubMedCfg(tool="chembl-da-test", email="team@ebi.ac.uk")

    text, error = pubmed.fetch_pubmed_batch(session, ["123", "456"], 0.0, cfg=cfg)

    assert error == ""
    assert isinstance(text, str)
    assert session.calls, "Expected PubMed request to be issued"
    url, timeout, kwargs = session.calls[-1]
    assert url.endswith("/efetch.fcgi")
    params = kwargs.get("params")
    assert isinstance(params, dict)
    assert params["tool"] == "chembl-da-test"
    assert params["email"] == "team@ebi.ac.uk"
    assert params["id"] == "123,456"
    assert timeout == (cfg.timeout_connect, cfg.timeout_read)


@pytest.mark.unit
def test_pubmed_cfg__rejects_placeholder_email() -> None:
    with pytest.raises(ValueError):
        PubMedCfg(email="user@example.org")


@pytest.mark.unit
def test_pubmed_cfg__requires_non_empty_tool() -> None:
    with pytest.raises(ValueError):
        PubMedCfg(tool="   ")


@pytest.mark.unit
def test_do_request__connect_timeout_logs_without_traceback(caplog: pytest.LogCaptureFixture) -> None:
    class _TimeoutSession:
        def get(
            self,
            url: str,
            *,
            timeout: float | tuple[float, float],
            **kwargs: object,
        ) -> None:
            raise requests.exceptions.ConnectTimeout("boom")

    session = _TimeoutSession()

    with caplog.at_level("INFO", logger="chembl"):
        data, error = pubmed._do_request(
            session,
            "https://api.openalex.org/works/pmid:20143779",
            0.0,
            retries=0,
            timeout=1.0,
        )

    assert data is None
    assert error == "boom"
    assert caplog.records, "Expected log record for failed request"
    record = caplog.records[-1]
    assert record.levelname == "WARNING"
    assert "request_fail" in record.getMessage()
    assert record.exc_info is None


@pytest.mark.unit
def test_do_request__deterministic_retry_delays(monkeypatch: pytest.MonkeyPatch) -> None:
    responses = [
        _SequencedResponse(503, text="Service unavailable"),
        _SequencedResponse(503, text="Try again"),
        _SequencedResponse(200, text="payload"),
    ]
    session = _SequencedSession(responses)
    recorded_delays: list[float] = []

    def _capture_sleep(value: float) -> None:
        recorded_delays.append(value)

    monkeypatch.setattr(pubmed, "sleep", _capture_sleep)

    retry_cfg = RetryCfg(backoff_factor=0.5, backoff_cap=None, jitter_seed=11)
    base_delay = 0.25

    data, error = pubmed._do_request(
        session,
        "https://example.org/resource",
        base_delay,
        expect_json=False,
        retries=2,
        timeout=1.0,
        retry_cfg=retry_cfg,
    )

    assert error == ""
    assert data == "payload"
    assert len(recorded_delays) == 2

    jitter = retry_cfg.build_jitter()
    assert jitter is not None
    expected_delays = [
        pubmed._retry_delay(1, base_delay, retry_cfg, timeout=1.0, jitter=jitter),
        pubmed._retry_delay(2, base_delay, retry_cfg, timeout=1.0, jitter=jitter),
    ]

    assert recorded_delays == pytest.approx(expected_delays)
def test_session_with_retry__disables_urllib3_retries() -> None:
    api_cfg = ApiCfg()
    retry_cfg = RetryCfg(max_attempts=4, backoff_factor=1.0, backoff_cap=None)

    with session_with_retry(api_cfg, retry_cfg) as session:
        adapter = session.get_adapter("https://")
        max_retries = adapter.max_retries

    assert isinstance(max_retries, urllib3_retry.Retry)
    assert max_retries.total == 0
    assert max_retries.connect == 0
    assert max_retries.read == 0
    assert max_retries.redirect == 0
    assert max_retries.status == 0
