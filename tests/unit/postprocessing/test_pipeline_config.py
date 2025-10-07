"""Unit tests for declarative post-processing pipeline configuration."""

from __future__ import annotations

import pytest

from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)


def test_load_pipeline_config__resolves_step_callables(monkeypatch):
    """The configuration loader returns callables defined in the target module."""

    monkeypatch.delenv("CHEMBL_ACTIVITY_PIPELINE_VERSION", raising=False)
    cfg = load_pipeline_config("activities")

    step_names = [step.name for step in cfg.steps]
    assert step_names == [
        "normalize_activity_records",
        "enrich_activity_quality",
        "finalize_activity_records",
    ]

    from library.postprocess.activities import steps as activity_steps

    assert cfg.steps[0].definition.func is activity_steps.normalize_activity_records
    assert cfg.steps[1].definition.func is activity_steps.enrich_activity_quality
    assert cfg.steps[2].definition.func is activity_steps.finalize_activity_records
    assert normalize_pipeline_version(cfg.pipeline_version) is None
    assert cfg.params["defaults"]["log_level"] == "INFO"


def test_load_pipeline_config__applies_env_overrides(monkeypatch):
    """Environment markers inside the YAML are expanded with defaults."""

    monkeypatch.setenv("CHEMBL_DOCUMENT_PIPELINE_VERSION", "2024.1")
    monkeypatch.setenv("POSTPROCESS_LOG_LEVEL", "DEBUG")

    cfg = load_pipeline_config("documents")

    assert cfg.pipeline_version == "2024.1"
    assert cfg.params["defaults"]["log_level"] == "DEBUG"


@pytest.mark.parametrize(
    "raw, expected",
    [
        (None, None),
        ("", None),
        ("auto", None),
        ("default", None),
        ("  2024.1  ", "2024.1"),
    ],
)
def test_normalize_pipeline_version__cases(raw: str | None, expected: str | None) -> None:
    """``normalize_pipeline_version`` trims sentinels and whitespace."""

    assert normalize_pipeline_version(raw) == expected
