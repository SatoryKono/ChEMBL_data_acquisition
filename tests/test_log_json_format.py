"""Tests for the shared logger instance in :mod:`library.log`."""

from __future__ import annotations

import io
import json
import sys

from library import log
from library.cli import LoggerConfig, configure_logger


def test_logger_emits_required_fields() -> None:
    """Configured logger produces JSON with ``ts``, ``level`` and ``event``."""

    buffer = io.StringIO()
    configure_logger(LoggerConfig(level="INFO", run_id="rid", stream=buffer))
    log.logger.info("test_event")
    record = json.loads(buffer.getvalue().splitlines()[0])
    configure_logger(LoggerConfig(stream=sys.stdout))
    assert record["event"] == "test_event"
    assert record["level"] == "INFO"
    assert "ts" in record
