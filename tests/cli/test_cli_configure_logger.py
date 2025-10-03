"""Tests for :func:`library.cli.configure_logger`."""

from __future__ import annotations

import io
import json

import pytest

from library.cli import LoggerConfig, configure_logger


def _parse(output: str) -> list[dict[str, object]]:
    return [json.loads(line) for line in output.strip().splitlines() if line]


def test_configure_logger_emits_json_without_overrides() -> None:
    """Default invocation should emit structured JSON logs."""

    stream = io.StringIO()
    logger = configure_logger(LoggerConfig(level="INFO", run_id="rid", stream=stream))
    logger.info("cli_event", foo="bar")
    records = _parse(stream.getvalue())
    assert len(records) == 1
    record = records[0]
    assert record["event"] == "cli_event"
    assert record["run_id"] == "rid"
    assert record["foo"] == "bar"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"fmt": "text"},
        {"datefmt": "%H:%M"},
        {"fmt": "text", "datefmt": "%H:%M"},
    ],
)
def test_configure_logger_rejects_format_override(kwargs: dict[str, str]) -> None:
    """Passing ``fmt`` or ``datefmt`` should raise ``ValueError``."""

    stream = io.StringIO()
    with pytest.raises(ValueError):
        configure_logger(LoggerConfig(stream=stream), **kwargs)
