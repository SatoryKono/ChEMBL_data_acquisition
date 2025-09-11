"""Tests that warnings are redirected to the structured logger.

ruff tests/test_warning_capture.py
black tests/test_warning_capture.py
mypy tests/test_warning_capture.py
pytest tests/test_warning_capture.py
"""

from __future__ import annotations

import json
import subprocess
import sys


def _parse(out: str) -> list[dict[str, object]]:
    """Return JSON objects parsed from newline-delimited ``out``."""
    return [json.loads(line) for line in out.strip().splitlines() if line]


def test_warnings_are_logged() -> None:
    """``warnings.warn`` emits a structured log record."""

    code = (
        "import sys, warnings, io; "
        "from library.logging_setup import LoggerConfig, configure_logger; "
        "buf = io.StringIO(); "
        "configure_logger(LoggerConfig(level='WARNING', run_id='rid', stream=buf)); "
        "warnings.warn('problem occurred'); "
        "print(buf.getvalue())"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    record = _parse(result.stdout)[0]
    assert record["level"] == "WARNING"
    assert record["event"] == "py.warnings"
    assert "problem occurred" in str(record.get("msg"))
    assert record["run_id"] == "rid"
