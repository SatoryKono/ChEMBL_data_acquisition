"""Regex-based tests for log timestamp format."""

from __future__ import annotations

import io
import json
import re
import sys

from library import log
from library.cli import LoggerConfig, configure_logger


def test_timestamp_has_iso_format() -> None:
    """Logger emits ISO 8601 timestamps with timezone."""
    buffer = io.StringIO()
    configure_logger(LoggerConfig(level="INFO", run_id="rid", stream=buffer))
    log.logger.info("sample_event")
    record = json.loads(buffer.getvalue().splitlines()[0])
    configure_logger(LoggerConfig(stream=sys.stdout))
    assert re.match(
        r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)?\+00:00",
        record["ts"],
    )
