from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_help_executes() -> None:
    """The script should be executable directly via Python."""
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "library.utils.cli_tools.get_input_initialisation",
            "--help",
        ],
        capture_output=True,
        text=True,
        check=False,
        cwd=Path(__file__).resolve().parents[1],
    )
    assert result.returncode == 0, result.stderr
    assert "Path to same document workbook" in result.stdout


def test_help_executes_file() -> None:
    """The script should run when invoked via its file path."""
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "library" / "utils" / "cli_tools" / "get_input_initialisation.py"
    result = subprocess.run(
        [sys.executable, str(script), "--help"],
        capture_output=True,
        text=True,
        check=False,
        cwd=repo_root,
    )
    assert result.returncode == 0, result.stderr
    assert "Path to same document workbook" in result.stdout
