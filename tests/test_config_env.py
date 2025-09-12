from __future__ import annotations

import sys
from pathlib import Path

from library.config import load_config
from library.logging_setup import LoggerConfig, configure_logger


def test_env_overrides_yaml_and_is_redacted(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text("api:\n  api_key: yaml-secret\n", encoding="utf-8")
    env_path = tmp_path / ".env"
    env_path.write_text("CHEMBL_DA__API__API_KEY=env-secret\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    cfg = load_config(config_path)
    assert cfg.api.api_key == "env-secret"

    test_logger = configure_logger(LoggerConfig(stream=sys.stdout))
    test_logger.info("test", api_key=cfg.api.api_key)
    output = capsys.readouterr().out
    assert "env-secret" not in output
    assert "***" in output
    monkeypatch.delenv("CHEMBL_DA__API__API_KEY", raising=False)
