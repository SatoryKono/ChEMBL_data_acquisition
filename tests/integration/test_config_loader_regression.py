from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import pytest

import library.config.loader as config_loader
from library.config import load_config


@pytest.mark.integration
def test_load_config__default_configuration_roundtrip(monkeypatch, tmp_path):
    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(tmp_path))
    monkeypatch.setattr(config_loader, "configure_rate_limiters", lambda _cfg: None)

    cfg, metadata = load_config(include_metadata=True)

    assert cfg.sources.chembl.api.chembl_base == "https://www.ebi.ac.uk/chembl/api/data"
    assert cfg.sources.chembl.api.rps == 20
    assert cfg.io.output_dir == (tmp_path / "output").resolve()

    rps_metadata = metadata.get("sources.chembl.api.rps")
    assert rps_metadata["value"] == 20
    assert rps_metadata["source"] == "config"
    assert rps_metadata["detail"] == config_loader.DEFAULT_CONFIG_RELATIVE.as_posix()
    assert Path(rps_metadata["detail"]).name == "config.yaml"


@pytest.mark.integration
def test_load_config__local_override_precedence(tmp_path, monkeypatch):
    base_cfg = tmp_path / "chembl.yaml"
    base_cfg.write_text(
        dedent(
            """
            local:
              io:
                output_dir: base_output
                cache_dir: base_cache
            system:
              log:
                level: INFO
            """
        ).strip()
    )

    local_cfg = base_cfg.with_name("chembl.local.yaml")
    local_cfg.write_text(
        dedent(
            """
            local:
              io:
                output_dir: override_output
            system:
              log:
                level: DEBUG
            """
        ).strip()
    )

    monkeypatch.setattr(config_loader, "configure_rate_limiters", lambda _cfg: None)

    cfg, metadata = load_config(base_cfg, include_metadata=True)

    assert cfg.io.output_dir == (base_cfg.parent / "override_output").resolve()
    assert cfg.io.cache_dir == (base_cfg.parent / "base_cache").resolve()
    assert cfg.log.level == "DEBUG"

    log_level_meta = metadata.get("system.log.level")
    assert log_level_meta["value"] == "DEBUG"
    assert log_level_meta["source"] == "config"
    assert Path(log_level_meta["detail"]).resolve() == local_cfg.resolve()
