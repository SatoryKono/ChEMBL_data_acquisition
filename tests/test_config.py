"""Tests for configuration loading and validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from library.config import Config, ConfigError, load_config


DATA_DIR = Path(__file__).resolve().parent / "data"


def test_load_config_applies_values(tmp_path: Path) -> None:
    """Values from the YAML file should populate the dataclass."""

    cfg_file = tmp_path / "config.yaml"
    cfg_file.write_text(
        f"""
api:
  chembl_base_url: https://example.com/api
timeouts:
  connect: 5
rate_limits:
  max_requests_per_second: 10
output:
  data_dir: {tmp_path / 'out'}
"""
    )

    cfg = load_config(cfg_file)
    assert cfg.api.chembl_base_url == "https://example.com/api"
    assert cfg.timeouts.connect == 5
    assert cfg.rate_limits.max_requests_per_second == 10
    assert cfg.output.data_dir == tmp_path / "out"


def test_load_config_missing_returns_defaults(tmp_path: Path) -> None:
    """Missing config files should return default settings."""

    cfg = load_config(tmp_path / "missing.yaml")
    assert isinstance(cfg, Config)
    assert cfg.timeouts.connect == 10.0


def test_invalid_url_raises() -> None:
    """Invalid URLs in the YAML file should raise ``ConfigError``."""

    with pytest.raises(ConfigError):
        load_config(DATA_DIR / "invalid_url.yaml")


def test_negative_timeout_raises() -> None:
    """Negative timeout values should raise ``ConfigError``."""

    with pytest.raises(ConfigError):
        load_config(DATA_DIR / "negative_timeout.yaml")


def test_environment_overrides(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Environment variables should override YAML values."""

    cfg_file = tmp_path / "config.yaml"
    cfg_file.write_text("timeouts:\n  connect: 5\n")

    monkeypatch.setenv("CHEMBL_TIMEOUTS__CONNECT", "15")
    monkeypatch.setenv("CHEMBL_OUTPUT__DATA_DIR", str(tmp_path / "env_out"))

    cfg = load_config(cfg_file)
    assert cfg.timeouts.connect == 15
    assert cfg.output.data_dir == tmp_path / "env_out"


def test_non_writable_directory_raises(tmp_path: Path) -> None:
    """Directories specified as files should raise ``ConfigError``."""

    bad_file = tmp_path / "file"
    bad_file.write_text("x")

    cfg_file = tmp_path / "config.yaml"
    cfg_file.write_text(f"output:\n  data_dir: {bad_file}\n")

    with pytest.raises(ConfigError):
        load_config(cfg_file)
