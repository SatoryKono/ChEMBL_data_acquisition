"""Tests for structured logging and retry utilities."""

from __future__ import annotations

import logging
from collections.abc import Callable
from typing import Any

import pytest
import requests

from library.clients.chembl_client import ChemblClient
from library.clients.crossref_client import CrossrefClient
from library.clients.pubchem_client import PubChemClient
from library.clients.uniprot_client import UniProtClient
from library.utils import logging as logging_utils
from library.utils.retry import with_retry


class _FakeResponse:
    def __init__(self, payload: dict[str, object], status_code: int = 200) -> None:
        self._payload = payload
        self.status_code = status_code

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.exceptions.HTTPError(response=None)

    def json(self) -> dict[str, object]:
        return self._payload


class _FlakySession:
    def __init__(self, failures: int) -> None:
        self._remaining_failures = failures
        self.calls: list[dict[str, object | None]] = []

    def get(
        self,
        url: str,
        *,
        params: dict[str, object] | None = None,
        timeout: float | None = None,
    ) -> _FakeResponse:
        self.calls.append({"url": url, "params": params, "timeout": timeout})
        if self._remaining_failures > 0:
            self._remaining_failures -= 1
            raise requests.exceptions.ConnectionError("boom")
        return _FakeResponse({"url": url, "params": params, "timeout": timeout})


def test_structured_logging_includes_context(monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture) -> None:
    logger = logging_utils.get_logger("tests.logging")
    times = iter([10.0, 10.5])
    monkeypatch.setattr(logging_utils, "_now", lambda: next(times))

    with caplog.at_level(logging.INFO):
        with logging_utils.log_context(run_id="run-1", stage="download", request_id="R42"):
            logger.info("download_complete", status="ok")

    record = caplog.records[0]
    data = getattr(record, "structured_data", {})
    assert data["run_id"] == "run-1"
    assert data["stage"] == "download"
    assert data["request_id"] == "R42"
    assert data["status"] == "ok"
    assert data["duration_s"] == 0.5
    formatter = logging_utils.StructuredFormatter()
    formatted = formatter.format(record)
    assert (
        formatted
        == "[INFO] [tests.logging] download_complete duration_s=0.5 request_id=R42 run_id=run-1 stage=download status=ok"
    )


def test_retry_decorator_retries_until_success(caplog: pytest.LogCaptureFixture) -> None:
    attempts = {"count": 0}
    logger = logging_utils.get_logger("tests.retry").bind(operation="demo")

    @with_retry(max_tries=3, timeout=0.1, logger=logger, log_event="demo_call")
    def flaky(*, timeout: float) -> str:
        attempts["count"] += 1
        assert timeout == 0.1
        if attempts["count"] < 3:
            raise requests.exceptions.ConnectionError("boom")
        return "ok"

    with caplog.at_level(logging.WARNING):
        result = flaky()

    assert result == "ok"
    assert attempts["count"] == 3
    text = caplog.text
    assert text.count("demo_call_retry") == 2


def test_retry_decorator_logs_fatal_error(caplog: pytest.LogCaptureFixture) -> None:
    attempts = {"count": 0}
    logger = logging_utils.get_logger("tests.retry.fail")

    @with_retry(max_tries=2, timeout=0.2, logger=logger, log_event="fatal_call")
    def always_fail(*, timeout: float) -> None:
        attempts["count"] += 1
        assert timeout == 0.2
        raise requests.exceptions.Timeout("timeout")

    with pytest.raises(requests.exceptions.Timeout):
        with caplog.at_level(logging.WARNING):
            always_fail()

    assert attempts["count"] == 2
    message = caplog.text
    assert "fatal_call_retry" in message
    assert "fatal_call_giveup" in message


@pytest.mark.parametrize(
    ("factory", "identifier", "expected_suffix"),
    (
        (ChemblClient, "CHEMBL1", "molecule/CHEMBL1"),
        (PubChemClient, "123", "compound/cid/123/JSON"),
        (UniProtClient, "P12345", "P12345.json"),
        (CrossrefClient, "10.1000/xyz", "10.1000/xyz"),
    ),
)
def test_clients_retry_and_timeout(
    caplog: pytest.LogCaptureFixture,
    factory: Callable[..., Any],
    identifier: str,
    expected_suffix: str,
) -> None:
    session = _FlakySession(failures=1)
    client = factory(session=session, timeout=0.3, max_tries=3)
    method_name = {
        ChemblClient: "get_molecule",
        PubChemClient: "get_compound",
        UniProtClient: "get_entry",
        CrossrefClient: "get_metadata",
    }[type(client)]
    method = getattr(client, method_name)

    with caplog.at_level(logging.INFO):
        payload = method(identifier)

    assert isinstance(payload, dict)
    assert len(session.calls) == 2
    last_call = session.calls[-1]
    assert last_call["timeout"] == 0.3
    assert expected_suffix in str(last_call["url"])
    assert caplog.text.count("_retry") == 1
