"""Smoke test for csv_utils_main runtime logging."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def test_csv_utils_main_logs_runtime(tmp_path: Path) -> None:
    """The CLI should log its execution duration."""
    input_csv = Path(__file__).parent / "data" / "csv_utils_input.csv"
    output_csv = tmp_path / "out.csv"
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "library.utils.cli_tools.csv_utils_main",
            "--input",
            str(input_csv),
            "--final-out",
            str(output_csv),
            "--key-cols",
            "a",
            "--log-level",
            "INFO",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    record = json.loads(proc.stdout.splitlines()[-1])
    assert record.get("event") == "run_completed"
    assert output_csv.exists()
