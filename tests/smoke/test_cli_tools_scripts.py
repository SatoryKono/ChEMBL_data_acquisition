"""Smoke tests for CLI modules under :mod:`library.utils.cli_tools`."""

from __future__ import annotations

import subprocess
import sys


def test_scripts_check_determinism_importable() -> None:
    """Running ``python -m library.utils.cli_tools.check_determinism`` should succeed."""

    subprocess.run(
        [
            sys.executable,
            "-m",
            "library.utils.cli_tools.check_determinism",
            "--log-level",
            "INFO",
        ],
        check=True,
    )
