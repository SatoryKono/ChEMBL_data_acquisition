"""Unit tests for declarative post-processing pipeline configuration."""

from __future__ import annotations

import importlib
from textwrap import dedent

import pytest

from library.postprocessing.pipeline.common.config import (
    PipelineConfigError,
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

    from library.postprocessing.pipeline.activities import steps as activity_steps

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
    "config_name,module_path,env_var,expected_steps,extra_sections",
    [
        (
            "activities",
            "library.postprocessing.pipeline.activities.steps",
            "CHEMBL_ACTIVITY_PIPELINE_VERSION",
            (
                (
                    "normalize_activity_records",
                    {
                        "relation_normalization": True,
                        "enforce_uppercase_units": True,
                    },
                ),
                (
                    "enrich_activity_quality",
                    {
                        "quality_terms": ["valid", "expert curated"],
                        "default_quality_flag": False,
                    },
                ),
                (
                    "finalize_activity_records",
                    {
                        "enforce_schema": True,
                        "numeric_identifier_dtype": "Int64",
                    },
                ),
            ),
            ("quality",),
        ),
        (
            "assays",
            "library.postprocessing.pipeline.assays.steps",
            "CHEMBL_ASSAY_PIPELINE_VERSION",
            (
                (
                    "normalize_assay_metadata",
                    {
                        "uppercase_categories": True,
                        "strip_whitespace": True,
                    },
                ),
                (
                    "enrich_assay_flags",
                    {
                        "confirmatory_terms": ["confirm", "primary"],
                        "default_flag": False,
                    },
                ),
                (
                    "finalize_assay_records",
                    {
                        "enforce_schema": True,
                        "normalize_identifiers": True,
                    },
                ),
            ),
            ("flags",),
        ),
        (
            "targets",
            "library.postprocessing.pipeline.targets.steps",
            "CHEMBL_TARGET_PIPELINE_VERSION",
            (
                (
                    "normalize_target_fields",
                    {
                        "normalize_taxonomy": True,
                        "fill_missing_identifiers": True,
                    },
                ),
                (
                    "enrich_target_synonyms",
                    {
                        "synonym_sources": ["chembl", "gtopdb"],
                        "preferred_separator": "; ",
                    },
                ),
                (
                    "finalize_target_records",
                    {
                        "enforce_schema": True,
                        "sort_by": ["target_chembl_id"],
                    },
                ),
            ),
            ("enrichment",),
        ),
        (
            "documents",
            "library.postprocessing.pipeline.documents.steps",
            "CHEMBL_DOCUMENT_PIPELINE_VERSION",
            (
                (
                    "normalize_document_fields",
                    {
                        "trim_whitespace": True,
                        "normalise_unicode": True,
                    },
                ),
                (
                    "enrich_document_publication_year",
                    {
                        "fallback_year": 1900,
                        "prefer_doi_year": True,
                    },
                ),
                (
                    "finalize_document_records",
                    {
                        "enforce_schema": True,
                        "ensure_unique_ids": True,
                    },
                ),
            ),
            ("enrichment",),
        ),
    ],
)
def test_pipeline_configs__resolve_steps_and_defaults(
    monkeypatch,
    config_name,
    module_path,
    env_var,
    expected_steps,
    extra_sections,
) -> None:
    """Each domain config resolves callables and exposes deterministic defaults."""

    monkeypatch.delenv(env_var, raising=False)
    monkeypatch.delenv("POSTPROCESS_LOG_LEVEL", raising=False)
    monkeypatch.delenv("POSTPROCESS_DEFAULT_ENCODING", raising=False)
    monkeypatch.delenv("POSTPROCESS_DEFAULT_CSV_SEPARATOR", raising=False)

    module = importlib.import_module(module_path)
    cfg = load_pipeline_config(config_name)

    assert [step.name for step in cfg.steps] == [name for name, _ in expected_steps]
    for step, (expected_name, expected_params) in zip(
        cfg.steps, expected_steps, strict=False
    ):
        assert step.name == expected_name
        assert step.definition.func is getattr(module, expected_name)
        assert step.params == expected_params

    assert cfg.params["defaults"]["log_level"] == "INFO"
    assert cfg.params["io"]["encoding"] == "utf-8"
    assert cfg.params["io"]["csv_sep"] == ","
    for section in ("defaults", "io", *extra_sections):
        assert section in cfg.params


def test_load_pipeline_config__errors_on_unsupported_step_parameters(tmp_path):
    """A clear error is raised when a step lists unsupported parameters."""

    config_text = dedent(
        """
        pipeline_version: "auto"
        enabled_steps:
          - name: normalize_document_fields
            callable: "library.postprocessing.pipeline.documents.steps:normalize_document_fields"
            params:
              unexpected: true
        """
    )
    config_path = tmp_path / "documents_invalid.yaml"
    config_path.write_text(config_text, encoding="utf-8")

    with pytest.raises(
        PipelineConfigError,
        match="step 'normalize_document_fields' defines unsupported parameters: unexpected",
    ):
        load_pipeline_config("documents", path=config_path)


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
def test_normalize_pipeline_version__cases(
    raw: str | None, expected: str | None
) -> None:
    """``normalize_pipeline_version`` trims sentinels and whitespace."""

    assert normalize_pipeline_version(raw) == expected
