"""Unit tests for :mod:`scripts.cleanup_project`."""

from __future__ import annotations

import os
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from scripts import cleanup_project


@pytest.mark.unit
def test_cleanup_project_check_detects_stale_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    repo_root = tmp_path / "repo"
    cache_dir = repo_root / "cache"
    cache_dir.mkdir(parents=True)
    stale_file = cache_dir / "stale.lock"
    stale_file.write_text("lock", encoding="utf-8")

    old_timestamp = (datetime.now(timezone.utc) - timedelta(days=5)).timestamp()
    os.utime(stale_file, (old_timestamp, old_timestamp))

    monkeypatch.setattr(cleanup_project, "PROJECT_ROOT", repo_root)
    monkeypatch.setattr(cleanup_project, "_load_tracked_paths", lambda _root: set())

    exit_code = cleanup_project.main(
        ["--check", "--categories", "pipeline-cache", "--days", "1"]
    )

    assert exit_code == 1
    assert stale_file.exists()
    captured = capsys.readouterr().out
    assert "Planned removals (1):" in captured
    assert "cache/stale.lock" in captured


@pytest.mark.unit
def test_cleanup_project_respects_age_threshold(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    repo_root = tmp_path / "repo"
    tmp_dir = repo_root / "temp"
    tmp_dir.mkdir(parents=True)
    recent_file = tmp_dir / "artifact.tmp"
    recent_file.write_text("temp", encoding="utf-8")

    monkeypatch.setattr(cleanup_project, "PROJECT_ROOT", repo_root)
    monkeypatch.setattr(cleanup_project, "_load_tracked_paths", lambda _root: set())

    exit_code = cleanup_project.main(["--categories", "tmp", "--days", "3"])

    assert exit_code == 0
    assert recent_file.exists()
    captured = capsys.readouterr().out
    assert "modified less than 3 day(s) ago" in captured
