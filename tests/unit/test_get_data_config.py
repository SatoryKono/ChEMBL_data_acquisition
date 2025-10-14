"""Tests for configuration overrides in ``scripts/get_data.py``."""

from __future__ import annotations

import re
from pathlib import Path
from types import SimpleNamespace

import pytest

from library.cli.commands import get_data
from scripts import get_data as compat_get_data


@pytest.mark.unit
def test_testitem_config_mapping__targets_nested_pipeline() -> None:
    """Ensure CLI overrides touch the nested testitem pipeline section."""

    expected_prefix = "sources.chembl.pipelines.testitem."
    dotted_targets = [
        ".".join(target) if isinstance(target, tuple) else str(target)
        for target in get_data._TESTITEM_CONFIG_MAPPING.values()
    ]
    assert all(target.startswith(expected_prefix) for target in dotted_targets)


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


@pytest.mark.unit
def test_stage_order_matches_readme_documentation() -> None:
    """The orchestrator stage order must match the README contract."""

    project_root = Path(__file__).resolve().parents[2]
    readme_path = project_root / "README.md"
    readme_text = readme_path.read_text(encoding="utf-8")

    pattern = re.compile(
        r"(Document)\s*→\s*(Target)\s*→\s*(Assay)\s*→\s*(Test item)\s*→\s*(Activity)"
    )
    match = pattern.search(readme_text)
    assert match is not None, "expected README to document the orchestrator stage order"

    documented_sequence = [value.lower().replace(" ", "") for value in match.groups()]
    stage_names = [stage.name for stage in compat_get_data.STAGES]

    assert stage_names == documented_sequence
