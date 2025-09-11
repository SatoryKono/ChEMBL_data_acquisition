from __future__ import annotations

import subprocess
import sys


def test_determinism_script_returns_zero() -> None:
    """Run determinism check script and assert a zero exit code.

    The script writes a small DataFrame twice using the deterministic CSV writer
    and compares their SHA-256 hashes. When the output is deterministic the
    script exits with code ``0``.
    """

    result = subprocess.run(
        [sys.executable, "scripts/check_determinism.py", "--log-level", "DEBUG"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
