"""Shared exports for the test item pipeline packages."""

from __future__ import annotations

from importlib import import_module

from .enrichment import enrich

_CATALOG_EXPORTS: tuple[str, ...] = (
    "LoadMoleculeHierarchyLookup",
    "PARENT_LOOKUP_SOURCE_CACHE",
    "PARENT_LOOKUP_SOURCE_LOOKUP",
    "PARENT_LOOKUP_SOURCE_PARTIAL",
    "PARENT_LOOKUP_SOURCE_SKIPPED",
    "PARENT_LOOKUP_SOURCE_SYNC",
    "ParentEnrichmentPreparation",
    "ParentEnrichmentResult",
    "ParentLookupPreparedData",
    "ParentLookupStats",
    "_DEFAULT_CATALOG_CFG",
    "_MOLECULE_HIERARCHY_COLUMNS",
    "_TYPO_PARENT_COLUMN",
    "_merge_parent_stats",
    "attach_parent_molecule_ids",
    "ensure_no_parant_column",
    "load_molecule_hierarchy_lookup",
    "load_parent_catalog",
    "molecule_catalog",
    "prepare_parent_enrichment",
    "query_parent_catalog",
    "run_parent_enrichment",
    "update_parent_catalog_cache",
    "write_parent_catalog_cache",
)

_CLI_EXPORTS: tuple[str, ...] = (
    "ReadInputIdsResult",
    "StageExecutionBudget",
    "StageWatchdog",
    "TestitemFetchError",
    "TestitemPipelineOptions",
    "TestitemPipelineStageError",
    "_FETCH_ERROR_SAMPLE_SIZE",
    "_log_missing_identifier_summary",
    "_prepare_pubchem_api_cfg",
    "apply_testitem_enrichment",
    "fetch_testitems",
    "finalize_output",
    "integrate_missing_identifiers",
    "read_input_ids",
    "run_testitem_pipeline",
)

_PUBCHEM_EXPORTS: tuple[str, ...] = (
    "PUBCHEM_CID_CACHE_ENCODING",
    "PUBCHEM_COLUMNS",
    "_CID_CACHE_MISSING",
    "_PUBCHEM_CACHE_SCHEMA_VERSION",
    "_PUBCHEM_SESSION_LOCK",
    "_PUBCHEM_SESSION_SIGNATURE",
    "_load_pubchem_cid_cache",
    "_merge_pubchem_properties",
    "_prepare_pubchem_caches",
    "_prefetch_parents",
    "_pubchem_session_signature",
    "_resolve_pubchem_cids",
    "_write_pubchem_cid_cache",
    "add_pubchem_data",
    "augment_pubchem",
    "resolve_pubchem_cid",
)

__all__ = [
    *_CATALOG_EXPORTS,
    *_CLI_EXPORTS,
    *_PUBCHEM_EXPORTS,
    "enrich",
]

_catalog_module = import_module(".catalog", __package__)
_cli_module = import_module(".cli", __package__)
_pubchem_module = import_module(".pubchem", __package__)

for _export_name in _CATALOG_EXPORTS:
    globals()[_export_name] = getattr(_catalog_module, _export_name)

for _export_name in _CLI_EXPORTS:
    globals()[_export_name] = getattr(_cli_module, _export_name)

for _export_name in _PUBCHEM_EXPORTS:
    globals()[_export_name] = getattr(_pubchem_module, _export_name)

del _catalog_module, _cli_module, _pubchem_module, _export_name
