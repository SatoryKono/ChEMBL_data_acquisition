
"""Structured logging tests for :mod:`library.logging_setup`.

Example
-------
Run linters and tests on this module with::

    ruff tests/test_logging_setup.py
    black tests/test_logging_setup.py
    mypy tests/test_logging_setup.py
    pytest tests/test_logging_setup.py
"""

from __future__ import annotations

import json
import sys
from typing import Any

import pytest

from library.logging_setup import LoggerConfig, configure_logger


def _parse(out: str) -> list[dict[str, Any]]:
    """Return JSON objects parsed from ``out`` separated by newlines."""

    return [json.loads(line) for line in out.strip().splitlines() if line]


def test_log_emits_json_per_line(capfd: pytest.CaptureFixture[str]) -> None:
    """Each call to :meth:`Logger.log` yields one JSON line with required fields."""

    logger = configure_logger(
        LoggerConfig(level="INFO", run_id="run123", stream=sys.stdout)
    )
    logger.log("INFO", "first")
    logger.log("INFO", "second")
    lines = _parse(capfd.readouterr().out)
    assert len(lines) == 2
    for rec in lines:
        assert rec["run_id"] == "run123"
        assert rec["level"] == "INFO"
        assert rec["event"] in {"first", "second"}
        assert "ts" in rec


def test_secret_redaction(capfd: pytest.CaptureFixture[str]) -> None:
    """Keys ending with sensitive suffixes are redacted."""

    logger = configure_logger(
        LoggerConfig(level="INFO", run_id="rid", stream=sys.stdout)
    )
    logger.info("secret", api_token="abc", db_password="xyz")
    record = _parse(capfd.readouterr().out)[0]
    assert record["api_token"] == "***"
    assert record["db_password"] == "***"


def test_stage_context_manager_success(
    capfd: pytest.CaptureFixture[str],
) -> None:
    """Successful stage logs start and done with matching context and elapsed."""

    logger = configure_logger(
        LoggerConfig(level="INFO", run_id="rid", stream=sys.stdout)
    )
    with logger.stage("extract"):
        pass
    start, done = _parse(capfd.readouterr().out)
    assert start["event"] == "extract_start"
    assert done["event"] == "extract_done"
    assert start["stage"] == done["stage"] == "extract"
    assert start["run_id"] == done["run_id"] == "rid"
    assert "elapsed" in done and isinstance(done["elapsed"], float)


def test_stage_context_manager_failure(
    capfd: pytest.CaptureFixture[str],
) -> None:
    """Exceptions inside a stage log a failure event."""

    logger = configure_logger(
        LoggerConfig(level="INFO", run_id="rid", stream=sys.stdout)
    )
    with pytest.raises(RuntimeError):
        with logger.stage("load"):
            raise RuntimeError("boom")
    start, fail = _parse(capfd.readouterr().out)
    assert start["event"] == "load_start"
    assert fail["event"] == "load_fail"
    assert start["stage"] == fail["stage"] == "load"
    assert start["run_id"] == fail["run_id"] == "rid"
