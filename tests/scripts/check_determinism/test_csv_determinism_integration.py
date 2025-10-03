from __future__ import annotations

import json
import subprocess
import sys


def test_csv_determinism_integration() -> None:
    """Run the check_determinism script and verify hashes match."""
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "library.utils.cli_tools.check_determinism",
            "--log-level",
            "DEBUG",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    records = [json.loads(line) for line in proc.stdout.splitlines()]
    hashes = [r["value"] for r in records if r.get("event") == "hash"]
    assert len(hashes) == 3
    assert len(set(hashes)) == 1
