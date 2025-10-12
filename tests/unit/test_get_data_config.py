"""Tests for configuration overrides in ``scripts/get_data.py``."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from library.cli.commands import get_data


@pytest.mark.unit
def test_testitem_config_mapping__targets_nested_pipeline() -> None:
    """Ensure CLI overrides touch the nested testitem pipeline section."""

    expected_prefix = "sources.chembl.pipelines.testitem."
    assert all(
        target.startswith(expected_prefix)
        for target in get_data._TESTITEM_CONFIG_MAPPING.values()
    )


@pytest.mark.unit
def test_load_pipeline_config__applies_testitem_limit_override(tmp_path: Path) -> None:
    """Loading the pipeline config should apply CLI overrides without errors."""

    project_root = Path(__file__).resolve().parents[2]
    cfg = SimpleNamespace(
        base_path=tmp_path,
        config_path=project_root / "config" / "config.yaml",
        limit=5,
        timeout=None,
        column=None,
        batch_size=None,
        offset=None,
    )

    config = get_data._load_pipeline_config(cfg, cfg.config_path)

    assert config.sources.chembl.pipelines.testitem.limit == 5
