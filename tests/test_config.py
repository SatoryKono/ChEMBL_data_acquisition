"""Tests for configuration loading and validation."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from library.config import (
    APISettings,
    Config,
    ConfigError,
    OutputPaths,
    TimeoutSettings,
    load_config,
)


def test_alias_expansion(tmp_path: Path) -> None:
    """Short-form aliases in YAML map to canonical keys."""

    cfg_file = tmp_path / "cfg.yaml"
    cfg_file.write_text(
        textwrap.dedent(
            """
            api:
              chembl:
                url: https://example.org/api
            """
        )
    )
    cfg = load_config(cfg_file)
    assert cfg.api.chembl_base_url == "https://example.org/api"


def test_env_var_override(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Environment variables override YAML values."""

    cfg_file = tmp_path / "cfg.yaml"
    cfg_file.write_text("timeouts:\n  connect: 5")
    monkeypatch.setenv("TIMEOUT_CONNECT", "99")
    cfg = load_config(cfg_file)
    assert cfg.timeouts.connect == 99.0


def test_missing_file_returns_defaults(tmp_path: Path) -> None:
    """Loading a missing file yields default settings."""

    path = tmp_path / "missing.yaml"
    cfg = load_config(path)
    assert cfg.api.chembl_base_url == APISettings().chembl_base_url


def test_negative_timeout_raises() -> None:
    """Validation rejects negative timeout values."""

    with pytest.raises(ConfigError):
        Config(timeouts=TimeoutSettings(connect=-1, read=10))


def test_invalid_url_raises() -> None:
    """Validation fails on malformed URLs."""

    api = APISettings(chembl_base_url="not-a-url")
    with pytest.raises(ConfigError):
        Config(api=api)


def test_non_writable_directory_raises(tmp_path: Path) -> None:
    """Non-writable output directories are rejected."""

    bad_path = tmp_path / "file"
    bad_path.write_text("x")
    with pytest.raises(ConfigError):
        Config(output=OutputPaths(data_dir=bad_path))
