from __future__ import annotations

import pytest

import library.config.loader as loader
from pathlib import Path

from library.config import Config, load_config


@pytest.mark.unit
def test_resolve_config_path__relative_name_uses_config_dir(tmp_path, monkeypatch):
    cfg_dir = tmp_path / "config"
    cfg_dir.mkdir()
    config_path = cfg_dir / "config.yaml"
    config_path.write_text("foo: bar\n", encoding="utf-8")

    monkeypatch.setattr(loader, "CONFIG_DIR", cfg_dir)
    monkeypatch.setattr(loader, "DEFAULT_CONFIG_PATH", config_path)
    monkeypatch.setattr(loader, "_DEFAULT_CONFIG_NAME", config_path.name)

    resolved = loader.resolve_config_path(config_path.name)

    assert resolved == config_path.resolve()


@pytest.mark.unit
def test_load_yaml_config__rejects_non_mapping(tmp_path, monkeypatch):
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text("- value\n- another\n", encoding="utf-8")

    monkeypatch.setattr(loader, "CONFIG_DIR", tmp_path)
    monkeypatch.setattr(loader, "DEFAULT_CONFIG_PATH", cfg_path)
    monkeypatch.setattr(loader, "_DEFAULT_CONFIG_NAME", cfg_path.name)

    with pytest.raises(loader.ConfigLoaderError):
        loader.load_yaml_config(cfg_path)


@pytest.mark.unit
def test_load_config__resolves_relative_paths_and_calls_runtime(tmp_path, monkeypatch):
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        """
        local:
          io:
            output_dir: output
            cache_dir: cache
        """.strip()
    )

    called: dict[str, Config] = {}

    def fake_configure(cfg):
        called["cfg"] = cfg

    monkeypatch.setattr(loader, "configure_rate_limiters", fake_configure)
    monkeypatch.setattr(
        "library.config.models.get_resource_path", lambda name: tmp_path / name
    )
    monkeypatch.setattr(
        "library.config.models.resolve_resource_reference",
        lambda value: Path(value),
    )

    cfg = load_config(cfg_path)

    assert cfg.io.output_dir == (cfg_path.parent / "output").resolve()
    assert cfg.io.cache_dir == (cfg_path.parent / "cache").resolve()
    assert called["cfg"] is cfg


@pytest.mark.unit
def test_load_config__metadata_records_cli_override(tmp_path, monkeypatch):
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text("system:\n  log:\n    level: INFO\n")

    monkeypatch.setattr(loader, "configure_rate_limiters", lambda cfg: None)
    monkeypatch.setattr(
        "library.config.models.get_resource_path", lambda name: tmp_path / name
    )
    monkeypatch.setattr(
        "library.config.models.resolve_resource_reference",
        lambda value: Path(value),
    )

    cfg, metadata = load_config(
        cfg_path,
        cli_overrides={"system.log.level": "DEBUG"},
        include_metadata=True,
        cli_sources={("system", "log", "level"): "--log-level"},
    )

    assert cfg.system.log.level == "DEBUG"
    cli_entry = metadata.cli_entry("--log-level")
    assert cli_entry == {"value": "DEBUG", "source": "cli", "detail": "--log-level"}
    path_entry = metadata.get("system.log.level")
    assert path_entry["source"] == "cli"
