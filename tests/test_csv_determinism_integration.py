from __future__ import annotations

import re
import subprocess
import sys


def test_csv_determinism_integration() -> None:
    """Run the check_determinism script and verify hashes match."""
    proc = subprocess.run(
        [sys.executable, "-m", "scripts.check_determinism", "--log-level", "DEBUG"],
        capture_output=True,
        text=True,
        check=True,
    )
    match1 = re.search(r"First hash: ([0-9a-f]{64})", proc.stderr)
    match2 = re.search(r"Second hash: ([0-9a-f]{64})", proc.stderr)
    assert match1 and match2
    assert match1.group(1) == match2.group(1)
