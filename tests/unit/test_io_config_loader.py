from __future__ import annotations

from pathlib import Path

import pytest

from library.io.config_loader import Config, load_config


@pytest.fixture(name="config_file")
def _config_file(tmp_path: Path) -> Path:
    cfg = tmp_path / "config.yaml"
    cfg.write_text("section:\n  key: 1\n", encoding="utf-8")
    return cfg


def test_load_config_returns_dataclass(config_file: Path) -> None:
    config = load_config(config_file)
    assert isinstance(config, Config)
    assert config.path == config_file
    assert config.to_dict()["section"]["key"] == 1


def test_environment_override_applied(config_file: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CHEMBL_DA__SECTION__KEY", "42")
    config = load_config(config_file)
    section = config.section("section")
    assert section["key"] == 42
