"""Smoke tests for thin CLI proxy modules in :mod:`scripts`."""

from __future__ import annotations

import subprocess
import sys


def test_scripts_check_determinism_importable() -> None:
    """Running ``python -m scripts.check_determinism`` should succeed."""

    subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.check_determinism",
            "--log-level",
            "INFO",
        ],
        check=True,
    )
