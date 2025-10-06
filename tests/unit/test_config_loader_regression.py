"""Regression tests for configuration loading."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from library.config.loader import load_config
from library.config.models import Config


@pytest.fixture
def clear_config_env(monkeypatch: pytest.MonkeyPatch) -> None:
    for key in list(os.environ):
        if key.startswith("CHEMBL_DA"):
            monkeypatch.delenv(key, raising=False)


@pytest.fixture
def minimal_config_path(
    tmp_path: Path, clear_config_env: None, monkeypatch: pytest.MonkeyPatch
) -> Path:
    def _fake_resource(name: str, base_dir: Path | None = None) -> Path:
        return tmp_path / name

    def _fake_resolve(value: Path | str) -> Path:
        return Path(value)

    monkeypatch.setattr("library.config.models.get_resource_path", _fake_resource)
    monkeypatch.setattr(
        "library.config.models.resolve_resource_reference", _fake_resolve
    )

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
        sources:
          chembl:
            molecule_catalog:
              hierarchy_lookup_path: null
        local:
          resources:
            dictionary_dir: "/tmp/dictionaries"
            iuphar_target_csv: "/tmp/iuphar_target.csv"
            iuphar_family_csv: "/tmp/iuphar_family.csv"
            uniprot_data_dir: "/tmp/uniprot"
            targets_type_csv: "/tmp/targets_type.csv"
        system:
          log:
            level: INFO
        """,
        encoding="utf-8",
    )
    return config_path


def test_load_config__defaults(minimal_config_path: Path) -> None:
    cfg = load_config(minimal_config_path)
    assert isinstance(cfg, Config)
    assert cfg.sources.chembl.api.timeout_read == 60.0
    assert cfg.system.log.level.upper() == "INFO"


def test_load_config__metadata_for_cli_overrides(minimal_config_path: Path) -> None:
    cli_path = ("system", "log", "level")
    cfg, metadata = load_config(
        minimal_config_path,
        cli_overrides={"system.log.level": "DEBUG"},
        cli_sources={cli_path: "--log-level"},
        include_metadata=True,
    )
    assert cfg.system.log.level == "DEBUG"
    snapshot = metadata.get("system.log.level")
    assert snapshot == {"value": "DEBUG", "source": "cli", "detail": "--log-level"}
