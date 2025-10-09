"""Unit tests for the CLI logging helpers."""

from __future__ import annotations

from pathlib import Path

from library.cli.logging import setup_cli_logging
from library.common.logging_setup import LoggerConfig


def test_setup_cli_logging__normalises_script_name(tmp_path: Path) -> None:
    """Log file prefixes should be derived from a normalised script name."""

    log_cfg = LoggerConfig(level="INFO", run_id="test-run")
    with setup_cli_logging("foo/bar.py", log_cfg, date_str="20240102", log_dir=tmp_path) as ctx:
        assert ctx.log_path.exists()
        assert ctx.log_path.name == "bar_20240102.log"

    second_cfg = LoggerConfig(level="INFO", run_id="second")
    with setup_cli_logging(" messy script name ", second_cfg, date_str="20240102", log_dir=tmp_path) as ctx:
        assert ctx.log_path.exists()
        assert ctx.log_path.name == "messy_script_name_20240102.log"
