"""Tests for :mod:`library.pipelines.configuration`."""

from __future__ import annotations

from library.pipelines.common.metadata import get_pipeline_version
from library.pipelines.configuration import (
    list_pipeline_configs,
    load_pipeline_config,
)


def test_list_pipeline_configs_contains_known_tables() -> None:
    names = set(list_pipeline_configs())
    expected = {
        "activity",
        "assay",
        "cellline",
        "document",
        "target",
        "testitem",
        "tissue",
    }
    assert expected.issubset(names)


def test_activity_pipeline_config_parameters_and_steps() -> None:
    cfg = load_pipeline_config("activity")
    assert cfg.pipeline_version.column == "pipeline_version"
    assert cfg.pipeline_version.value == get_pipeline_version()
    mapping = cfg.parameters.mapping_for()
    assert mapping["timeout"] == "activity.timeout"
    assert mapping["dry_run"] == "activity.dry_run"
    step_names = {step.name for step in cfg.steps}
    assert "fetch_activities" in step_names
    assert "finalize_export" in step_names


def test_document_pipeline_modes_and_shared_parameters() -> None:
    cfg = load_pipeline_config("document")
    shared = cfg.parameters.mapping_for()
    assert shared["openalex_rps"] == "openalex.rps"
    pubmed_mapping = cfg.parameters.mapping_for(mode="pubmed")
    assert pubmed_mapping["sleep"] == "document.pubmed.sleep"
    all_mapping = cfg.parameters.mapping_for(mode="all")
    assert all_mapping["chembl_chunk_size"] == "document.all.chunk_size"
    step_lookup = {step.name: step for step in cfg.steps}
    assert "enrich_literature_partners" in step_lookup
    assert "pubmed" in step_lookup["enrich_literature_partners"].applies_to


def test_target_pipeline_commands_mapping() -> None:
    cfg = load_pipeline_config("target")
    chembl_map = cfg.parameters.mapping_for(command="chembl")
    assert chembl_map["chunk_size"] == "target.chembl.chunk_size"
    all_map = cfg.parameters.mapping_for(command="all")
    assert all_map["chembl_chunk_size"] == "target.chembl.chunk_size"
    assert all_map["target_csv"] == "target.all.target_csv"
