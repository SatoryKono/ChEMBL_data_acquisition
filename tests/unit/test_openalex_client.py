from __future__ import annotations

from urllib.parse import quote

import pytest
import requests
import responses

from library.clients.openalex import fetch_openalex
from library.config.models import OpenAlexCfg

pytestmark = pytest.mark.unit


class _DummyLimiter:
    def __init__(self) -> None:
        self.calls: int = 0

    def acquire(self) -> None:
        self.calls += 1


@pytest.fixture(autouse=True)
def disable_network() -> None:  # pragma: no cover - override global fixture
    """Allow mocked HTTP requests via ``responses``."""

    yield


@pytest.fixture()
def cfg() -> OpenAlexCfg:
    return OpenAlexCfg(mailto="chembl-testing@ebi.ac.uk")


@pytest.fixture()
def limiter() -> _DummyLimiter:
    return _DummyLimiter()


def _expected_url(cfg: OpenAlexCfg, pmid: str) -> str:
    base = cfg.base.rstrip("/")
    return f"{base}/works/pmid:{pmid}?mailto={quote(cfg.mailto)}"


def test_fetch_openalex__success(cfg: OpenAlexCfg, limiter: _DummyLimiter) -> None:
    pmid = "20143779"
    url = _expected_url(cfg, pmid)

    with requests.Session() as session:
        with responses.RequestsMock(assert_all_requests_are_fired=True) as rs:
            rs.add(rs.GET, url, json={"id": "https://openalex.org/W20143779"}, status=200)

            data, error = fetch_openalex(session, pmid, cfg=cfg, limiter=limiter)

    assert limiter.calls == 1
    assert error == ""
    assert isinstance(data, dict)
    assert data["id"] == "https://openalex.org/W20143779"


@pytest.mark.parametrize(
    "status, body, expectation",
    [
        (404, {"error": "Not found"}, "PMID not found"),
        (200, "not json", "Invalid JSON:"),
    ],
)
def test_fetch_openalex__errors_return_none(
    status: int,
    body: object,
    expectation: str,
    cfg: OpenAlexCfg,
    limiter: _DummyLimiter,
    caplog: pytest.LogCaptureFixture,
) -> None:
    pmid = "18737329"
    url = _expected_url(cfg, pmid)

    with requests.Session() as session:
        with responses.RequestsMock(assert_all_requests_are_fired=True) as rs:
            if isinstance(body, str):
                rs.add(rs.GET, url, body=body, status=status, content_type="application/json")
            else:
                rs.add(rs.GET, url, json=body, status=status)

            with caplog.at_level("INFO", logger="chembl"):
                data, error = fetch_openalex(session, pmid, cfg=cfg, limiter=limiter)

    assert limiter.calls == 1
    assert data is None
    assert error.startswith(expectation)
    assert any("request_fail" in record.getMessage() for record in caplog.records)


def test_fetch_openalex__timeout_logged(
    cfg: OpenAlexCfg,
    limiter: _DummyLimiter,
    caplog: pytest.LogCaptureFixture,
) -> None:
    pmid = "18737330"
    url = _expected_url(cfg, pmid)
    cfg = cfg.copy(update={"retries": 0})

    with requests.Session() as session:
        with responses.RequestsMock(assert_all_requests_are_fired=True) as rs:
            rs.add(rs.GET, url, body=requests.exceptions.Timeout("simulated timeout"))

            with caplog.at_level("WARNING", logger="chembl"):
                data, error = fetch_openalex(session, pmid, cfg=cfg, limiter=limiter)

    assert limiter.calls == 1
    assert data is None
    assert error == "simulated timeout"
    assert any(
        record.levelname == "WARNING" and "request_fail" in record.getMessage()
        for record in caplog.records
    )
