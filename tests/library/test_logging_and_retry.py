"""Tests covering structured logging and retry helpers."""

from __future__ import annotations

import io
import logging
from dataclasses import dataclass
from typing import Any

import pytest
import requests

from library.clients.chembl_client import ChemblServiceClient
from library.utils.logging import get_logger, log_context
from library.utils.retry import retryable


@pytest.mark.unit
def test_get_logger_produces_structured_message(capsys: pytest.CaptureFixture[str]) -> None:
    logger = get_logger("tests.logging")
    with log_context(run_id="run-1", stage="load"):
        logger.info("stage_complete", result="ok")

    output = capsys.readouterr().out.strip()
    assert output.startswith("[INFO] [tests.logging] stage_complete")
    assert "run_id='run-1'" in output
    assert "stage='load'" in output
    assert "duration_s=" in output
    assert "result='ok'" in output


class _RetryProbe:
    def __init__(self) -> None:
        self.calls = 0

    @retryable(max_attempts=3, timeout=2.5, logger=get_logger("tests.retry"))
    def execute(self, *, timeout: float) -> str:
        self.calls += 1
        assert timeout == 2.5
        if self.calls < 3:
            raise requests.exceptions.ConnectionError("boom")
        return "ok"


@pytest.mark.unit
def test_retryable_retries_until_success() -> None:
    buffer = io.StringIO()
    logger = get_logger("tests.retry")
    for handler in logger.handlers:
        logger.removeHandler(handler)
    stream = logging.StreamHandler(buffer)
    stream.setFormatter(logging.Formatter("[%(levelname)s] [%(name)s] %(message)s"))
    logger.addHandler(stream)

    probe = _RetryProbe()
    result = probe.execute()

    assert result == "ok"
    assert probe.calls == 3
    assert buffer.getvalue().count("http_retry") == 2


@pytest.mark.unit
def test_retryable_logs_giveup() -> None:
    buffer = io.StringIO()
    logger = get_logger("tests.retry.fail")
    for handler in logger.handlers:
        logger.removeHandler(handler)
    stream = logging.StreamHandler(buffer)
    stream.setFormatter(logging.Formatter("[%(levelname)s] [%(name)s] %(message)s"))
    logger.addHandler(stream)

    attempts: int = 0

    @retryable(max_attempts=2, logger=get_logger("tests.retry.fail"))
    def _always_fail(*, timeout: float) -> None:
        nonlocal attempts
        attempts += 1
        raise requests.exceptions.Timeout("timeout")

    with pytest.raises(requests.exceptions.Timeout):
        _always_fail()

    assert attempts == 2
    log_output = buffer.getvalue()
    assert "http_fail" in log_output


@dataclass
class _StubResponse:
    url: str
    status_code: int
    payload: dict[str, Any]

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.exceptions.HTTPError(
                f"{self.status_code} error for {self.url}", response=self
            )

    def json(self) -> dict[str, Any]:
        return self.payload


class _FlakySession:
    def __init__(self) -> None:
        self.calls: list[tuple[str, float | None]] = []
        self._step = 0

    def get(
        self, url: str, *, params: Any | None = None, timeout: float | None = None
    ) -> _StubResponse:
        del params
        self.calls.append((url, timeout))
        if self._step == 0:
            self._step += 1
            raise requests.exceptions.Timeout("temporary")
        return _StubResponse(url, 200, {"ok": True})


@pytest.mark.unit
def test_chembl_client_logs_context() -> None:
    session = _FlakySession()
    logger = get_logger("tests.chembl_client")
    buffer = io.StringIO()
    for handler in logger.handlers:
        logger.removeHandler(handler)
    stream = logging.StreamHandler(buffer)
    stream.setFormatter(logging.Formatter("[%(levelname)s] [%(name)s] %(message)s"))
    logger.addHandler(stream)
    client = ChemblServiceClient(
        base_url="https://example.test/chembl/api/data",
        session=session,
        timeout=5.0,
        max_attempts=3,
        run_id="run-chembl",
        logger=logger,
    )

    data = client.fetch_assay("CHEMBL1")

    assert data == {"ok": True}
    assert session.calls == [
        ("https://example.test/chembl/api/data/assay/CHEMBL1", 5.0),
        ("https://example.test/chembl/api/data/assay/CHEMBL1", 5.0),
    ]

    captured = buffer.getvalue()
    assert "run_id='run-chembl'" in captured
    assert "http_success" in captured
    assert "http_retry" in captured
    assert "chembl_fetch" in captured
