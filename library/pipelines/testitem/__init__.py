"""Test item pipeline orchestration helpers.

Changelog
========
* Add compatibility re-exports for the modern test item pipeline API.
* Re-export public helpers from :mod:`library.testitem_pipeline` for legacy
  imports.
* Document :mod:`library.pipelines.testitem` as the canonical module while the
  legacy :mod:`library.testitem_pipeline` shim continues to proxy to this
  package.
"""

from __future__ import annotations

# ===== Modules =====
from pathlib import Path

from library.clients import pubchem as pc
from library.common.csv_utils import (
    write_csv_chunks_deterministic as write_csv_deterministic,
)
from library.common.log import logger
from library.config import Config
from library.integration import molecule_catalog, pubchem_library as pl
from library.integration.chembl_client import ChemblClient
from library.integration.molecule_catalog import (
    load_parent_catalog,
    query_parent_catalog,
    update_parent_catalog_cache,
    write_parent_catalog_cache,
)
from library.metadata import file_sha256, write_meta_yaml
from library.pipelines.assay.chembl_assay import TESTITEM_PUBCHEM_COLUMNS
from library.pipelines.common import PipelineRunResult
from library.table_quality import analyze_table_quality
from library.validation import validate_testitems

from . import enrichment as testitem_enrichment
from .catalog import (
    LoadMoleculeHierarchyLookup,
    PARENT_LOOKUP_SOURCE_CACHE,
    PARENT_LOOKUP_SOURCE_LOOKUP,
    PARENT_LOOKUP_SOURCE_PARTIAL,
    PARENT_LOOKUP_SOURCE_SKIPPED,
    PARENT_LOOKUP_SOURCE_SYNC,
    ParentEnrichmentPreparation,
    ParentEnrichmentResult,
    ParentLookupPreparedData,
    ParentLookupStats,
    _DEFAULT_CATALOG_CFG,
    _MOLECULE_HIERARCHY_COLUMNS,
    _TYPO_PARENT_COLUMN,
    _merge_parent_stats,
    attach_parent_molecule_ids,
    ensure_no_parant_column,
    load_molecule_hierarchy_lookup,
    prepare_parent_enrichment,
    run_parent_enrichment,
)
from .cli import (
    ReadInputIdsResult,
    TestitemFetchError,
    TestitemPipelineOptions,
    TestitemPipelineStageError,
    StageExecutionBudget,
    StageWatchdog,
    _FETCH_ERROR_SAMPLE_SIZE,
    _log_missing_identifier_summary,
    _prepare_pubchem_api_cfg,
    apply_testitem_enrichment,
    fetch_testitems,
    finalize_output,
    integrate_missing_identifiers,
    read_input_ids,
    run_testitem_pipeline,
)
from .enrichment import enrich
from .pubchem import (
    PUBCHEM_CID_CACHE_ENCODING,
    PUBCHEM_COLUMNS,
    _CID_CACHE_MISSING,
    _PUBCHEM_CACHE_SCHEMA_VERSION,
    _PUBCHEM_SESSION_LOCK,
    _PUBCHEM_SESSION_SIGNATURE,
    _load_pubchem_cid_cache,
    _merge_pubchem_properties,
    _prepare_pubchem_caches,
    _prefetch_parents,
    _pubchem_session_signature,
    _resolve_pubchem_cids,
    _write_pubchem_cid_cache,
    add_pubchem_data,
    augment_pubchem,
    resolve_pubchem_cid,
)

__all__ = [
    "ChemblClient",
    "LoadMoleculeHierarchyLookup",
    "PARENT_LOOKUP_SOURCE_CACHE",
    "PARENT_LOOKUP_SOURCE_LOOKUP",
    "PARENT_LOOKUP_SOURCE_PARTIAL",
    "PARENT_LOOKUP_SOURCE_SKIPPED",
    "PARENT_LOOKUP_SOURCE_SYNC",
    "PUBCHEM_CID_CACHE_ENCODING",
    "PUBCHEM_COLUMNS",
    "ParentEnrichmentPreparation",
    "ParentEnrichmentResult",
    "ParentLookupPreparedData",
    "ParentLookupStats",
    "ReadInputIdsResult",
    "StageExecutionBudget",
    "StageWatchdog",
    "TESTITEM_PUBCHEM_COLUMNS",
    "TestitemFetchError",
    "TestitemPipelineOptions",
    "TestitemPipelineStageError",
    "_CID_CACHE_MISSING",
    "_DEFAULT_CATALOG_CFG",
    "_FETCH_ERROR_SAMPLE_SIZE",
    "_MOLECULE_HIERARCHY_COLUMNS",
    "_PUBCHEM_CACHE_SCHEMA_VERSION",
    "_PUBCHEM_SESSION_LOCK",
    "_PUBCHEM_SESSION_SIGNATURE",
    "_TYPO_PARENT_COLUMN",
    "_log_missing_identifier_summary",
    "_load_pubchem_cid_cache",
    "_merge_pubchem_properties",
    "_merge_parent_stats",
    "_prepare_pubchem_api_cfg",
    "_prepare_pubchem_caches",
    "_prefetch_parents",
    "_pubchem_session_signature",
    "_resolve_pubchem_cids",
    "_write_pubchem_cid_cache",
    "add_pubchem_data",
    "analyze_table_quality",
    "apply_testitem_enrichment",
    "attach_parent_molecule_ids",
    "augment_pubchem",
    "enrich",
    "ensure_no_parant_column",
    "fetch_testitems",
    "file_sha256",
    "finalize_output",
    "integrate_missing_identifiers",
    "load_molecule_hierarchy_lookup",
    "load_parent_catalog",
    "logger",
    "molecule_catalog",
    "pc",
    "pl",
    "prepare_parent_enrichment",
    "query_parent_catalog",
    "read_input_ids",
    "resolve_pubchem_cid",
    "run_parent_enrichment",
    "run_pipeline",
    "run_testitem_pipeline",
    "update_parent_catalog_cache",
    "validate_testitems",
    "write_csv_deterministic",
    "write_meta_yaml",
    "write_parent_catalog_cache",
]

__all__.append("testitem_enrichment")


def run_pipeline(config: Config, options: TestitemPipelineOptions) -> PipelineRunResult:
    """Execute the test item pipeline and return a common result structure."""

    cfg = config.model_copy(deep=True)
    pipelines = cfg.sources.chembl.pipelines
    section = pipelines.testitem
    updates: dict[str, object] = {}
    if getattr(options, "limit", None) is not None:
        updates["limit"] = options.limit
    if getattr(options, "offset", None) is not None:
        updates["offset"] = options.offset
    if updates:
        pipelines.testitem = section.model_copy(update=updates)

    output_candidate = getattr(options, "output_csv", None)
    if output_candidate is not None:
        output_path = Path(output_candidate)
    else:
        output_path = Path(options.input_csv)

    exit_code = run_testitem_pipeline(cfg, options)
    reason = None if exit_code == 0 else "pipeline_failed"
    written = None if exit_code != 0 else True
    return PipelineRunResult(
        exit_code=exit_code,
        output_path=output_path,
        executed=True,
        reason=reason,
        written=written,
    )
