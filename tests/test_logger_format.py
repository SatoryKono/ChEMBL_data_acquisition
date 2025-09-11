from __future__ import annotations

import pytest

from library.cli import LoggerConfig, configure_logger
from library.log import logger


def test_shared_logger_format(capfd: pytest.CaptureFixture[str]) -> None:
    """Shared logger uses configured format."""
    cfg = LoggerConfig(run_id="test", level="INFO")
    configure_logger(cfg, fmt="%(levelname)s|%(message)s")
    logger.info("hello")
    captured = capfd.readouterr()
    assert captured.err.strip() == "INFO|hello"
