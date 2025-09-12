from __future__ import annotations

import shutil
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

import library.metadata as metadata
from library.metadata import Stats, _git_sha, write_meta_yaml


def test_write_meta_yaml_creates_file(tmp_path: Path) -> None:
    data_src = Path(__file__).parent / "data" / "meta_input.csv"
    csv_path = tmp_path / "output.csv"
    shutil.copy(data_src, csv_path)

    stats: Stats = {
        "rows_total": 2,
        "rows_kept": 2,
        "rows_dropped": 0,
        "output_sha256": "deadbeef",
    }

    meta_path = write_meta_yaml(
        csv_path=csv_path,
        command="unit-test",
        config_subset={"api_key": "secret", "password": "topsecret"},
        inputs={"source": "dummy"},
        stats=stats,
        schema="TestSchema",
    )

    assert meta_path.exists()
    with meta_path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)

    assert data["command"] == "unit-test"
    assert data["config"]["api_key"] == "***"
    assert data["config"]["password"] == "***"
    assert data["inputs"]["source"] == "dummy"
    assert data["stats"] == stats
    assert data["schema"] == "TestSchema"

    required_keys = {
        "generated_at",
        "git_sha",
        "python_version",
        "platform",
        "command",
        "config",
        "inputs",
        "stats",
        "schema",
    }
    assert required_keys <= data.keys()


def test_git_sha_env_var(monkeypatch: pytest.MonkeyPatch) -> None:
    """_git_sha uses the ``GIT_SHA`` environment variable when available."""

    monkeypatch.setenv("GIT_SHA", "envsha")
    with patch("library.metadata.subprocess.check_output") as mock:
        assert _git_sha() == "envsha"
        mock.assert_not_called()


def test_git_sha_missing_git_executable(monkeypatch: pytest.MonkeyPatch) -> None:
    """_git_sha returns UNKNOWN and warns when git is absent."""

    monkeypatch.setattr(shutil, "which", lambda _cmd: None)
    messages: list[str] = []
    monkeypatch.setattr(
        metadata.logger, "warning", lambda msg, *args: messages.append(msg % args)
    )

    assert metadata._git_sha() == "UNKNOWN"
    assert "Git executable not found" in messages[0]


def test_git_sha_missing_git_dir(monkeypatch: pytest.MonkeyPatch) -> None:
    """_git_sha returns UNKNOWN and warns when .git directory is missing."""

    repo_root = Path(metadata.__file__).resolve().parent.parent
    original_exists = metadata.Path.exists

    def mock_exists(self: Path) -> bool:  # type: ignore[override]
        if self == repo_root / ".git":
            return False
        return original_exists(self)

    monkeypatch.setattr(metadata.Path, "exists", mock_exists)
    messages: list[str] = []
    monkeypatch.setattr(
        metadata.logger, "warning", lambda msg, *args: messages.append(msg % args)
    )

    assert metadata._git_sha() == "UNKNOWN"
    assert f"No .git directory found at {repo_root}" in messages[0]
