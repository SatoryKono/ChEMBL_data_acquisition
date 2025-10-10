from __future__ import annotations

from pathlib import Path

import pytest

import library.config.loader as loader
from config.paths import DICTIONARY_DIR
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


@pytest.mark.unit
def test_load_config__metadata_includes_env_warnings(tmp_path, monkeypatch):
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text("system:\n  log:\n    level: INFO\n", encoding="utf-8")

    monkeypatch.setenv("CHEMBL_DA__SYSTEM__UNKNOWN", "42")
    monkeypatch.setattr(loader, "configure_rate_limiters", lambda cfg: None)
    monkeypatch.setattr(
        "library.config.models.get_resource_path", lambda name: tmp_path / name
    )
    monkeypatch.setattr(
        "library.config.models.resolve_resource_reference", lambda value: Path(value)
    )

    _, metadata = load_config(cfg_path, include_metadata=True)

    assert metadata.env_warnings == [
        "CHEMBL_DA__SYSTEM__UNKNOWN: unknown configuration path 'system.unknown'"
    ]


@pytest.mark.unit
def test_load_config__preserves_dictionary_resource_reference(tmp_path, monkeypatch):
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        """
        local:
          resources:
            dictionary_dir: dictionary_root
        """.strip()
    )

    monkeypatch.setattr(loader, "configure_rate_limiters", lambda cfg: None)
    monkeypatch.setattr(
        loader,
        "list_resource_names",
        lambda *, validate=True, base_dir=None: ("dictionary_root",),
    )
    loader._dictionary_resource_names.cache_clear()
    monkeypatch.setattr(
        "library.config.models.resolve_resource_reference",
        lambda value: DICTIONARY_DIR if value == "dictionary_root" else Path(value),
    )
    monkeypatch.setattr(
        "library.config.models.get_resource_path",
        lambda name: DICTIONARY_DIR if name == "dictionary_root" else Path(name),
    )

    cfg = load_config(cfg_path)

    assert cfg.resources.dictionary_dir == DICTIONARY_DIR.resolve()


@pytest.mark.unit
def test_dictionary_resource_names__uses_manifest_without_validation(monkeypatch):
    calls: list[bool] = []

    def fake_list_resource_names(*, validate: bool = True, base_dir=None):
        calls.append(validate)
        return ("foo", "bar")

    monkeypatch.setattr(loader, "list_resource_names", fake_list_resource_names)
    loader._dictionary_resource_names.cache_clear()

    result = loader._dictionary_resource_names()

    assert result == frozenset({"foo", "bar"})
    assert calls == [False]
