"""Test item pipeline orchestration helpers.

Changelog
========
* Consolidate the public API under ``library.pipelines.testitem`` while
  preserving the legacy ``library.testitem_pipeline`` compatibility layer.
"""

from __future__ import annotations

import requests

from library.clients import pubchem as pc
from library.common.csv_utils import write_csv_chunks_deterministic as write_csv_deterministic
from library.common.log import logger
from library.integration import pubchem_library as pl
from library.integration.chembl_client import ChemblClient
from library.qa.validation import validate_testitems

from . import catalog as catalog_module
from . import enrichment as testitem_enrichment
from . import cli as cli_module
from . import pubchem as pubchem_module
from .enrichment import enrich

LoadMoleculeHierarchyLookup = catalog_module.LoadMoleculeHierarchyLookup
PARENT_LOOKUP_SOURCE_CACHE = catalog_module.PARENT_LOOKUP_SOURCE_CACHE
PARENT_LOOKUP_SOURCE_LOOKUP = catalog_module.PARENT_LOOKUP_SOURCE_LOOKUP
PARENT_LOOKUP_SOURCE_PARTIAL = catalog_module.PARENT_LOOKUP_SOURCE_PARTIAL
PARENT_LOOKUP_SOURCE_SKIPPED = catalog_module.PARENT_LOOKUP_SOURCE_SKIPPED
PARENT_LOOKUP_SOURCE_SYNC = catalog_module.PARENT_LOOKUP_SOURCE_SYNC
ParentEnrichmentPreparation = catalog_module.ParentEnrichmentPreparation
ParentEnrichmentResult = catalog_module.ParentEnrichmentResult
ParentLookupPreparedData = catalog_module.ParentLookupPreparedData
ParentLookupStats = catalog_module.ParentLookupStats
_DEFAULT_CATALOG_CFG = catalog_module._DEFAULT_CATALOG_CFG
_MOLECULE_HIERARCHY_COLUMNS = catalog_module._MOLECULE_HIERARCHY_COLUMNS
_TYPO_PARENT_COLUMN = catalog_module._TYPO_PARENT_COLUMN
_merge_parent_stats = catalog_module._merge_parent_stats
attach_parent_molecule_ids = catalog_module.attach_parent_molecule_ids
ensure_no_parant_column = catalog_module.ensure_no_parant_column
load_molecule_hierarchy_lookup = catalog_module.load_molecule_hierarchy_lookup
load_parent_catalog = catalog_module.load_parent_catalog
molecule_catalog = catalog_module.molecule_catalog
prepare_parent_enrichment = catalog_module.prepare_parent_enrichment
query_parent_catalog = catalog_module.query_parent_catalog
run_parent_enrichment = catalog_module.run_parent_enrichment
update_parent_catalog_cache = catalog_module.update_parent_catalog_cache
write_parent_catalog_cache = catalog_module.write_parent_catalog_cache

ReadInputIdsResult = cli_module.ReadInputIdsResult
TestitemFetchError = cli_module.TestitemFetchError
TestitemPipelineOptions = cli_module.TestitemPipelineOptions
TestitemPipelineStageError = cli_module.TestitemPipelineStageError
_FETCH_ERROR_SAMPLE_SIZE = cli_module._FETCH_ERROR_SAMPLE_SIZE
_log_missing_identifier_summary = cli_module._log_missing_identifier_summary
_prepare_pubchem_api_cfg = cli_module._prepare_pubchem_api_cfg
analyze_table_quality = cli_module.analyze_table_quality
apply_testitem_enrichment = cli_module.apply_testitem_enrichment
fetch_testitems = cli_module.fetch_testitems
finalize_output = cli_module.finalize_output
file_sha256 = cli_module.file_sha256
integrate_missing_identifiers = cli_module.integrate_missing_identifiers
read_input_ids = cli_module.read_input_ids
run_testitem_pipeline = cli_module.run_testitem_pipeline
write_meta_yaml = cli_module.write_meta_yaml

PUBCHEM_CID_CACHE_ENCODING = pubchem_module.PUBCHEM_CID_CACHE_ENCODING
PUBCHEM_COLUMNS = list(pubchem_module.PUBCHEM_COLUMNS)
_CID_CACHE_MISSING = pubchem_module._CID_CACHE_MISSING
_PUBCHEM_CACHE_SCHEMA_VERSION = pubchem_module._PUBCHEM_CACHE_SCHEMA_VERSION
_PUBCHEM_SESSION_LOCK = pubchem_module._PUBCHEM_SESSION_LOCK
_PUBCHEM_SESSION_SIGNATURE = pubchem_module._PUBCHEM_SESSION_SIGNATURE
_load_pubchem_cid_cache = pubchem_module._load_pubchem_cid_cache
_merge_pubchem_properties = pubchem_module._merge_pubchem_properties
_prepare_pubchem_caches = pubchem_module._prepare_pubchem_caches
_prefetch_parents = pubchem_module._prefetch_parents
_pubchem_session_signature = pubchem_module._pubchem_session_signature
_resolve_pubchem_cids = pubchem_module._resolve_pubchem_cids
_write_pubchem_cid_cache = pubchem_module._write_pubchem_cid_cache
add_pubchem_data = pubchem_module.add_pubchem_data
augment_pubchem = pubchem_module.augment_pubchem
resolve_pubchem_cid = pubchem_module.resolve_pubchem_cid

_CATALOG_EXPORTS = list(catalog_module.__all__)
_CLI_EXPORTS = [
    "ReadInputIdsResult",
    "TestitemFetchError",
    "TestitemPipelineOptions",
    "TestitemPipelineStageError",
    "_FETCH_ERROR_SAMPLE_SIZE",
    "_log_missing_identifier_summary",
    "_prepare_pubchem_api_cfg",
    "analyze_table_quality",
    "apply_testitem_enrichment",
    "fetch_testitems",
    "finalize_output",
    "file_sha256",
    "integrate_missing_identifiers",
    "read_input_ids",
    "run_testitem_pipeline",
    "write_meta_yaml",
]
_PUBCHEM_EXPORTS = list(pubchem_module.__all__)

__all__ = [
    "enrich",
    "testitem_enrichment",
]
__all__ += _CATALOG_EXPORTS
__all__ += _CLI_EXPORTS
__all__ += _PUBCHEM_EXPORTS

# Maintain deterministic order without duplicates.
__all__ = list(dict.fromkeys(__all__))
