from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any, cast
from unittest.mock import patch

import pytest
import yaml

import library.utils.git as git_utils
from library.metadata import Stats, write_meta_yaml


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


def test_write_meta_yaml_preserves_existing_columns(tmp_path: Path) -> None:
    csv_path = tmp_path / "output.csv"
    csv_path.write_text("id\n1\n", encoding="utf-8")
    meta_path = csv_path.with_name(csv_path.name + ".meta.yaml")
    original_meta = {
        "columns": ["id"],
        "dtypes": {"id": "string"},
        "schema": "Initial",
    }
    meta_path.write_text(yaml.safe_dump(original_meta), encoding="utf-8")

    stats: Stats = {
        "rows_total": 1,
        "rows_kept": 1,
        "rows_dropped": 0,
        "output_sha256": "feedface",
    }

    write_meta_yaml(
        csv_path=csv_path,
        command="unit-test",
        config_subset={},
        inputs={},
        stats=stats,
        schema="Updated",
    )

    data = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert data["columns"] == ["id"]
    assert data["dtypes"] == {"id": "string"}
    assert data["schema"] == "Updated"
    assert data["stats"] == stats


def test_git_sha_env_var(monkeypatch: pytest.MonkeyPatch) -> None:
    """_git_sha uses the ``GIT_SHA`` environment variable when available."""

    git_utils._git_sha.cache_clear()
    monkeypatch.setenv("GIT_SHA", "envsha")
    with patch("library.utils.git.subprocess.run") as mock:
        assert git_utils._git_sha() == "envsha"
        mock.assert_not_called()


def test_git_sha_missing_git_executable(monkeypatch: pytest.MonkeyPatch) -> None:
    """_git_sha returns UNKNOWN and warns when git is absent."""

    git_utils._git_sha.cache_clear()
    monkeypatch.setattr(git_utils, "_read_head_sha", lambda *_: None)
    monkeypatch.setattr(shutil, "which", lambda _cmd: None)
    messages: list[str] = []
    logger = cast(Any, getattr(git_utils, "logger"))  # noqa: B009
    monkeypatch.setattr(
        logger,
        "warning",
        lambda msg, *args, **kwargs: messages.append(msg % args),
    )

    assert git_utils._git_sha() == "UNKNOWN"
    assert "git_executable_missing" in messages[0]


def test_git_sha_missing_git_dir(monkeypatch: pytest.MonkeyPatch) -> None:
    """_git_sha returns UNKNOWN and warns when .git directory is missing."""

    git_utils._git_sha.cache_clear()
    monkeypatch.setattr(git_utils, "_resolve_git_dir", lambda *_: None)
    messages: list[str] = []
    logger = cast(Any, getattr(git_utils, "logger"))  # noqa: B009
    monkeypatch.setattr(
        logger,
        "warning",
        lambda msg, *args, **kwargs: messages.append(msg % args),
    )

    assert git_utils._git_sha() == "UNKNOWN"
    assert "git_directory_missing" in messages[0]
