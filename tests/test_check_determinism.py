"""Regression tests for :mod:`scripts.check_determinism`."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from scripts import check_determinism


class _DummyCompletedProcess(subprocess.CompletedProcess[str]):
    """Helper providing a string-specialised :class:`CompletedProcess`."""

    def __init__(self) -> None:
        super().__init__(args=["python", "scripts/get_activity_data.py"], returncode=0)


def _render_csv(limit: int) -> str:
    """Return deterministic CSV content for ``limit`` identifiers."""

    rows = "\n".join(str(index) for index in range(limit))
    return f"activity_id\n{rows}\n"


def test_main_succeeds(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Check that the determinism probe exits successfully."""

    outputs: dict[str, str] = {}

    def _fake_run_activity(limit: int, destination: Path) -> subprocess.CompletedProcess[str]:
        assert limit == 3
        content = _render_csv(limit)
        destination.write_text(content)
        outputs[destination.name] = content
        return _DummyCompletedProcess()

    monkeypatch.setattr(check_determinism, "_run_activity", _fake_run_activity)
    monkeypatch.setattr(check_determinism.tempfile, "mkdtemp", lambda prefix: str(tmp_path))

    exit_code = check_determinism.main(["--limit", "3"])

    assert exit_code == 0
    assert outputs["first.csv"] == outputs["second.csv"]

