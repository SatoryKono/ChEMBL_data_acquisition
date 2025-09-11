from __future__ import annotations

import subprocess
import sys


def test_check_determinism_script() -> None:
    """Ensure the check_determinism script exits successfully."""
    result = subprocess.run(
        [sys.executable, "-m", "scripts.check_determinism"], check=False
    )
    assert result.returncode == 0
