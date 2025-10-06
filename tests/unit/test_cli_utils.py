"""Unit tests for :mod:`library.cli_utils`."""

from __future__ import annotations

import argparse
from unittest.mock import Mock

import pytest

from library import cli_utils
from library.cli import LoggerConfig


@pytest.mark.unit
def test_run_cli_command__logs_exc_info_on_value_error(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure ``exc_info`` is attached when configuration resolution fails."""

    parser = argparse.ArgumentParser()
    args = argparse.Namespace(config=object(), verbose=False, log_level=None)
    log_cfg = LoggerConfig(level="INFO", run_id="test-run", redact_secrets=False)

    logger = Mock()
    monkeypatch.setattr(cli_utils.cli, "configure_logger", lambda cfg: logger)

    exit_code = cli_utils.run_cli_command(
        args=args,
        parser=parser,
        log_cfg=log_cfg,
        mapping={},
        run=lambda *_: 0,
        logger=logger,
    )

    assert exit_code == 1
    error_call = logger.error.call_args
    assert error_call is not None
    exc_info = error_call.kwargs.get("exc_info")
    assert exc_info is not None
    assert isinstance(exc_info, ValueError)
    assert error_call.kwargs["error"] == "configuration path must be provided"
