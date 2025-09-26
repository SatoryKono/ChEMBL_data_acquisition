"""Ensure that variables from a `.env` file override the YAML configuration."""

from __future__ import annotations

import pytest

from library.config import ConfigError, load_config


def test_env_file_overrides(tmp_path, monkeypatch) -> None:
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        "sources:\n  chembl:\n    api:\n      rps: 5\n",
        encoding="utf8",
    )

    env_path = tmp_path / ".env"
    env_path.write_text("CHEMBL_DA__SOURCES__CHEMBL__API__RPS=7\n", encoding="utf8")

    for line in env_path.read_text(encoding="utf8").splitlines():
        key, value = line.split("=", 1)
        monkeypatch.setenv(key, value)
    monkeypatch.setenv(
        "CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT",
        "test-agent/1.0 (mailto:test@example.org)",
    )

    cfg = load_config(cfg_path)
    assert cfg.api.rps == 7


def test_env_override_invalid_value(tmp_path, monkeypatch) -> None:
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text("", encoding="utf8")

    monkeypatch.setenv("CHEMBL_DA_RPS", "0")

    with pytest.raises(ConfigError, match="CHEMBL_DA_RPS must be ≥1"):
        load_config(cfg_path)


def test_env_override_invalid_yaml(tmp_path, monkeypatch) -> None:
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text("", encoding="utf8")

    monkeypatch.setenv("CHEMBL_DA_RPS", "[1, 2")

    with pytest.raises(ConfigError, match="CHEMBL_DA_RPS could not be parsed"):
        load_config(cfg_path)
