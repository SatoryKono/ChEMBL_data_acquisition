from __future__ import annotations

from pathlib import Path

import pytest

from library.config import Config, ConfigMetadata, load_config


@pytest.fixture
def sample_config(tmp_path: Path) -> Path:
    config_text = """
sources:
  chembl:
    api:
      chembl_base: "https://chembl.example/api"
      user_agent: "chembl-da/1.0 (mailto:chembl-data@ebi.ac.uk)"
local:
  resources:
    dictionary_dir: "$CHEMBL_DA_BASE_PATH/resources"
system:
  log:
    level: INFO
"""
    path = tmp_path / "config.yaml"
    path.write_text(config_text)
    return path


@pytest.mark.unit
def test_load_config__metadata_tracks_sources(sample_config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("CHEMBL_DA_BASE", raising=False)
    base_path = tmp_path / "data"
    base_path.mkdir()

    cfg, metadata = load_config(sample_config, base_path=base_path, include_metadata=True)

    assert isinstance(cfg, Config)
    assert isinstance(metadata, ConfigMetadata)
    expected_dir = (base_path / "resources").resolve()
    assert cfg.local.resources.dictionary_dir == expected_dir

    entry = metadata.get(("local", "resources", "dictionary_dir"))
    assert entry["source"] == "config"
    assert Path(entry["value"]) == expected_dir


@pytest.mark.unit
def test_load_config__applies_environment_alias(sample_config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    override_url = "https://override.example/api"
    monkeypatch.setenv("CHEMBL_DA_BASE", override_url)

    cfg, metadata = load_config(sample_config, base_path=tmp_path, include_metadata=True)

    assert cfg.sources.chembl.api.chembl_base == override_url
    entry = metadata.get(("sources", "chembl", "api", "chembl_base"))
    assert entry["source"] == "env"
    assert entry["value"] == override_url


@pytest.mark.unit
def test_load_config__records_cli_override(sample_config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("CHEMBL_DA_BASE", raising=False)
    override_dir = (tmp_path / "out").resolve()
    cli_overrides = {"local.io.output_dir": str(override_dir)}
    cli_sources = {("local", "io", "output_dir"): "output_dir"}

    cfg, metadata = load_config(
        sample_config,
        base_path=tmp_path,
        include_metadata=True,
        cli_overrides=cli_overrides,
        cli_sources=cli_sources,
    )

    assert cfg.local.io.output_dir == override_dir
    entry = metadata.get(("local", "io", "output_dir"))
    assert entry["source"] == "cli"
    assert entry["detail"] == "output_dir"
    assert Path(entry["value"]) == override_dir
