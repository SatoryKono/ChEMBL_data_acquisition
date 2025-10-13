"""Integration tests for the determinism smoke test offline mode."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.integration
@pytest.mark.pipeline_scenario("idempotence")
def test_check_determinism__offline_mode(tmp_path: Path) -> None:
    _ = tmp_path
    repo_root = Path(__file__).resolve().parents[2]
    fixtures_dir = repo_root / "tests" / "resources" / "expected_get_data"

    result = subprocess.run(
        [
            sys.executable,
            "scripts/check_determinism.py",
            "--limit",
            "3",
            "--offline",
            "--fixtures-dir",
            str(fixtures_dir),
        ],
        capture_output=True,
        text=True,
        cwd=repo_root,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "Deterministic output confirmed" in result.stdout
