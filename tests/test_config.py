"""Tests for the configuration loader."""

from pathlib import Path

import pytest

from library.config import ConfigError, load_config


def test_load_config_success(tmp_path: Path) -> None:
    cfg_path = tmp_path / "cfg.yaml"
    cfg_path.write_text("key: value\n", encoding="utf8")
    cfg = load_config(cfg_path)
    assert cfg == {"key": "value"}


def test_load_config_missing(tmp_path: Path) -> None:
    missing = tmp_path / "missing.yaml"
    with pytest.raises(ConfigError) as excinfo:
        load_config(missing)
    assert "configuration file not found" in str(excinfo.value)


def test_load_config_invalid_yaml(tmp_path: Path) -> None:
    bad = tmp_path / "bad.yaml"
    bad.write_text("key: [1, 2", encoding="utf8")  # malformed YAML
    with pytest.raises(ConfigError) as excinfo:
        load_config(bad)
    assert "invalid YAML" in str(excinfo.value)
