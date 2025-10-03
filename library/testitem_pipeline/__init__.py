"""Test item pipeline helpers exposed as a public API."""

from __future__ import annotations

import json
import sys
from importlib import import_module
from importlib.util import (
    find_spec,
    module_from_spec,
    spec_from_file_location,
)
from pathlib import Path
from types import ModuleType
from typing import Any

import requests

_CATALOG_EXPORTS = (
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


def _load_local_module(module_name: str) -> ModuleType:
    """Load a sibling module directly from disk as a fallback."""

    package_name = __name__
    qualified_name = f"{package_name}.{module_name}"
    module_path = Path(__file__).with_name(f"{module_name}.py")
    spec = spec_from_file_location(qualified_name, module_path)
    if spec is None or spec.loader is None:
        raise ModuleNotFoundError(qualified_name)
    module = module_from_spec(spec)
    sys.modules[qualified_name] = module
    spec.loader.exec_module(module)
    return module


def _export_from_module(module: ModuleType, names: tuple[str, ...]) -> None:
    """Populate the module globals with the selected attributes."""

    for name in names:
        globals()[name] = getattr(module, name)
catalog_module_name = f"{__name__}.catalog"
catalog_spec = find_spec(catalog_module_name)

if catalog_spec is not None:
    catalog_module = import_module(catalog_module_name)
else:  # pragma: no cover - environment-specific fallback
    catalog_module = _load_local_module("catalog")

_export_from_module(catalog_module, _CATALOG_EXPORTS)


class _LazyModuleProxy:
    """Lazily import *module_name* on first attribute access."""

    __slots__ = ("_module", "_module_name")

    def __init__(self, module_name: str) -> None:
        self._module_name = module_name
        self._module: ModuleType | None = None

    def _ensure_loaded(self) -> ModuleType:
        if self._module is None:
            self._module = import_module(self._module_name)
        return self._module

    def __getattr__(self, name: str) -> Any:  # noqa: D401 - behave like the target module
        module = self._ensure_loaded()
        return getattr(module, name)

    def __dir__(self) -> list[str]:
        module = self._ensure_loaded()
        return sorted(set(dir(module)))


testitem_enrichment = _LazyModuleProxy("library.pipelines.testitem.enrichment")


from .cli import (
    ReadInputIdsResult,
    TestitemFetchError,
    TestitemPipelineOptions,
    TestitemPipelineStageError,
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
from library.integration import pubchem_library as pl
from library.integration.chembl_client import ChemblClient
from library.clients import pubchem as pc
from library.common.csv_utils import write_csv_chunks_deterministic as write_csv_deterministic
from library.common.log import logger
from library.metadata import file_sha256, write_meta_yaml
from library.table_quality import analyze_table_quality
from library.validation import validate_testitems

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
    "_load_pubchem_cid_cache",
    "_log_missing_identifier_summary",
    "_merge_parent_stats",
    "_merge_pubchem_properties",
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
    "ensure_no_parant_column",
    "fetch_testitems",
    "file_sha256",
    "finalize_output",
    "integrate_missing_identifiers",
    "json",
    "load_molecule_hierarchy_lookup",
    "load_parent_catalog",
    "logger",
    "molecule_catalog",
    "pc",
    "pl",
    "prepare_parent_enrichment",
    "query_parent_catalog",
    "read_input_ids",
    "requests",
    "resolve_pubchem_cid",
    "run_parent_enrichment",
    "run_testitem_pipeline",
    "testitem_enrichment",
    "update_parent_catalog_cache",
    "validate_testitems",
    "write_meta_yaml",
    "write_csv_deterministic",
    "write_parent_catalog_cache",
]
