"""Tests that warnings are redirected to the structured logger.

Example
-------
Run linters and tests on this module with::

    ruff tests/test_logging.py
    black tests/test_logging.py
    mypy tests/test_logging.py
    pytest tests/test_logging.py
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
    assert record["logger"] == "py.warnings"
    assert "problem occurred" in str(record.get("event"))
    assert record["run_id"] == "rid"
