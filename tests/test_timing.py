from __future__ import annotations

from io import StringIO
from time import perf_counter

from library.cli import LoggerConfig, configure_logger
from library.timing import log_duration


def test_log_duration_logs_and_returns_value() -> None:
    """Ensure ``log_duration`` logs the elapsed time in seconds."""
    stream = StringIO()
    configure_logger(LoggerConfig(stream=stream))
    start = perf_counter()
    duration = log_duration(start)
    assert duration >= 0
    assert "duration_sec" in stream.getvalue()
