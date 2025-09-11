from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_help_executes() -> None:
    """The script should be executable directly via Python."""
    result = subprocess.run(
        [sys.executable, "scripts/get_input_initialisation.py", "--help"],
        capture_output=True,
        text=True,
        check=False,
        cwd=Path(__file__).resolve().parents[1],
    )
    assert result.returncode == 0, result.stderr
    assert "Path to same document workbook" in result.stdout
