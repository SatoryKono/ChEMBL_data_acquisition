"""Unit tests for logging helpers used by CLI entrypoints."""

from __future__ import annotations

import pytest

from library.cli.logging import setup_cli_logging
from library.common.logging_setup import LoggerConfig


@pytest.mark.unit
def test_setup_cli_logging__normalises_script_name(tmp_path, monkeypatch):
    base_dir = tmp_path / "base"
    log_dir = base_dir / "logs"
    log_dir.mkdir(parents=True)

    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(base_dir))
    monkeypatch.setattr("library.cli.logging._DEFAULT_LOG_DIR", log_dir)
    monkeypatch.setattr("library.cli.logging._current_date_str", lambda: "20240102")

    cfg = LoggerConfig(level="INFO", run_id="run-42", handlers=[])

    script_identifier = "scripts\\get-activity data.py"

    with setup_cli_logging(script_identifier, cfg) as ctx:
        log_path = ctx.log_path
        assert log_path.parent == log_dir
        assert log_path.name == "get-activity_data_20240102.log"
        assert log_path.exists()

    # The helper should leave behind the log file for subsequent inspection.
    assert (log_dir / "get-activity_data_20240102.log").exists()


@pytest.mark.unit
def test_setup_cli_logging__keeps_existing_names(tmp_path, monkeypatch):
    base_dir = tmp_path / "workspace"
    log_dir = base_dir / "logs"

    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(base_dir))
    monkeypatch.setattr("library.cli.logging._DEFAULT_LOG_DIR", log_dir)
    monkeypatch.setattr("library.cli.logging._current_date_str", lambda: "20240102")

    cfg = LoggerConfig(level="INFO", run_id="run-42", handlers=[])

    with setup_cli_logging("get_activity_data", cfg) as ctx:
        assert ctx.log_path.name == "get_activity_data_20240102.log"

    assert (log_dir / "get_activity_data_20240102.log").exists()
