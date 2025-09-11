"""Ensure that variables from a `.env` file override the YAML configuration."""

from __future__ import annotations

from library.config import load_config


def test_env_file_overrides(tmp_path, monkeypatch) -> None:
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text("api:\n  rps: 5\n", encoding="utf8")

    env_path = tmp_path / ".env"
    env_path.write_text("CHEMBL_DA__API__RPS=7\n", encoding="utf8")

    for line in env_path.read_text(encoding="utf8").splitlines():
        key, value = line.split("=", 1)
        monkeypatch.setenv(key, value)
    monkeypatch.setenv(
        "CHEMBL_DA__API__USER_AGENT", "test-agent/1.0 (mailto:test@example.org)"
    )

    cfg = load_config(cfg_path)
    assert cfg.api.rps == 7
