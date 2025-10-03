from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from library.common.logging_setup import LoggerConfig, configure_logger
from library.utils.cli_tools import csv_utils_main


def test_cli_run_id(capfd: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    """CLI emits consistent run_id across run_start and run_done events."""
    input_csv = Path("tests/data/csv_utils_input.csv")
    output_csv = tmp_path / "out.csv"
    logger = configure_logger(LoggerConfig(run_id="test-run", stream=sys.stdout))
    with logger.stage("run"):
        exit_code = csv_utils_main.main(
            [
                "--input",
                str(input_csv),
                "--final-out",
                str(output_csv),
                "--key-cols",
                "a",
            ]
        )
    assert exit_code == 0
    captured = capfd.readouterr().out.splitlines()
    events = [json.loads(line) for line in captured if line.strip()]
    assert events[0]["event"] == "run_start"
    assert events[-1]["event"] == "run_done"
    run_id_start = events[0]["run_id"]
    run_id_done = events[-1]["run_id"]
    assert run_id_start == run_id_done
    assert run_id_start
