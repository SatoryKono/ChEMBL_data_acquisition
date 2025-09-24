from __future__ import annotations

from pathlib import Path

from scripts.check_determinism import run_check


def test_run_check_repeatable(tmp_path: Path) -> None:
    """Running the determinism check twice yields consistent results."""

    workdir = tmp_path / "check"
    workdir.mkdir()

    assert run_check(workdir)
    # Re-run in the same directory to ensure files are overwritten deterministically.
    assert run_check(workdir)
