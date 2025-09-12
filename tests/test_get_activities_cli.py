from __future__ import annotations

import subprocess
import sys


def test_get_activities_cli_smoke() -> None:
    """Run the get-activities CLI and expect a zero exit code."""
    result = subprocess.run(
        [sys.executable, "-m", "scripts.get_activities", "--limit", "1", "--dry-run"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
