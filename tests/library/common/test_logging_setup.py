"""Structured logging tests for :mod:`library.common.logging_setup`.

Example
-------
Run linters and tests on this module with::

    ruff tests/test_logging_setup.py
    black tests/test_logging_setup.py
    mypy tests/test_logging_setup.py
    pytest tests/test_logging_setup.py
"""

from __future__ import annotations

import io
import json
import logging
import sys
import threading
from typing import Any

import pytest

from library.common.logging_setup import LoggerConfig, configure_logger


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


def test_multithreaded_logging_emits_valid_json(
    capfd: pytest.CaptureFixture[str],
) -> None:
    """Concurrent logging should emit one JSON object per line."""

    logger = configure_logger(
        LoggerConfig(level="INFO", run_id="rid", stream=sys.stdout)
    )

    def worker(idx: int) -> None:
        logger.info("event", thread=idx)

    threads = [threading.Thread(target=worker, args=(i,)) for i in range(10)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    lines = capfd.readouterr().out.strip().splitlines()
    assert len(lines) == 10
    for line in lines:
        json.loads(line)


def test_logger_bind_preserves_parent_context(
    capfd: pytest.CaptureFixture[str],
) -> None:
    """Binding context should not mutate the original logger."""

    base = configure_logger(LoggerConfig(level="INFO", run_id="rid", stream=sys.stdout))
    child = base.bind(stage="extract")

    child.info("child_event")
    base.info("parent_event")

    records = _parse(capfd.readouterr().out)
    assert records[0]["event"] == "child_event"
    assert records[0]["stage"] == "extract"
    parent = next(rec for rec in records if rec["event"] == "parent_event")
    assert "stage" not in parent


def test_configure_logger_without_replacing_root_preserves_handlers() -> None:
    """Opting out of root replacement keeps existing handlers intact."""

    root = logging.getLogger()
    warnings_logger = logging.getLogger("py.warnings")
    structured_logger = logging.getLogger("library.structured")

    class DummyHandler(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:  # pragma: no cover - noop
            return None

    original_root_handlers = list(root.handlers)
    original_root_level = root.level
    original_warnings_handlers = list(warnings_logger.handlers)
    original_warnings_level = warnings_logger.level
    original_warnings_propagate = warnings_logger.propagate
    original_structured_handlers = list(structured_logger.handlers)
    original_structured_level = structured_logger.level
    original_structured_propagate = structured_logger.propagate

    root_handler = DummyHandler()
    warnings_handler = DummyHandler()
    root.handlers = [root_handler]
    warnings_logger.handlers = [warnings_handler]
    warnings_logger.propagate = True

    buffer = io.StringIO()

    try:
        logger = configure_logger(
            LoggerConfig(level="INFO", run_id="rid", stream=buffer),
            replace_root=False,
        )

        assert root.handlers == [root_handler]
        assert root.level == original_root_level
        assert warnings_logger.handlers == [warnings_handler]
        assert warnings_logger.level == original_warnings_level
        assert warnings_logger.propagate is True

        structured_logger.info("std_event")
        logger.info("direct_event")

        records = [json.loads(line) for line in buffer.getvalue().splitlines()]
        events = {record["event"] for record in records}
        assert {"std_event", "direct_event"}.issubset(events)
    finally:
        root.handlers = original_root_handlers
        root.setLevel(original_root_level)
        warnings_logger.handlers = original_warnings_handlers
        warnings_logger.setLevel(original_warnings_level)
        warnings_logger.propagate = original_warnings_propagate
        structured_logger.handlers = original_structured_handlers
        structured_logger.setLevel(original_structured_level)
        structured_logger.propagate = original_structured_propagate
