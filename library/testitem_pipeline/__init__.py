"""Test item pipeline helpers exposed as a public API."""

from __future__ import annotations

# ruff: noqa: E402
import itertools
import json
import sys
from importlib import import_module, resources
from importlib.machinery import SourcelessFileLoader
from importlib.util import (
    find_spec,
    module_from_spec,
    spec_from_file_location,
    spec_from_loader,
)
from pathlib import Path
from types import ModuleType
from typing import Any

import requests

# ===== changelog =====
# 2024-05-21: Enhance optional module loader to support compiled-only installs.
# 2024-05-22: Build exports list programmatically to avoid merge artifacts.
# 2025-03-09: Keep resource handles open while loading optional modules.

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

_CLI_EXPORTS = (
    "ReadInputIdsResult",
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

_PUBCHEM_EXPORTS = (
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


def _load_local_module(module_name: str) -> ModuleType:
    """Load a sibling module from package resources as a fallback."""

    package_name = __name__
    qualified_name = f"{package_name}.{module_name}"
    suffixes = (".py", ".pyc")

    def _candidate_roots() -> tuple[Path, ...]:
        """Return directories that may contain the requested module."""

        roots: list[Path] = []

        package = sys.modules.get(package_name)
        if package is not None:
            spec = getattr(package, "__spec__", None)
            if spec and getattr(spec, "submodule_search_locations", None):
                for location in spec.submodule_search_locations:
                    try:
                        roots.append(Path(location))
                    except TypeError:
                        continue

        roots.append(Path(__file__).resolve().parent)
        # Ensure deterministic order while removing duplicates.
        seen: set[Path] = set()
        ordered: list[Path] = []
        for root in roots:
            if root not in seen:
                seen.add(root)
                ordered.append(root)
        return tuple(ordered)

    def _load_from_path(module_path: Path) -> ModuleType:
        if module_path.suffix == ".pyc":
            loader = SourcelessFileLoader(qualified_name, str(module_path))
            spec = spec_from_loader(qualified_name, loader)
        else:
            spec = spec_from_file_location(qualified_name, str(module_path))
        if spec is None or spec.loader is None:
            raise ModuleNotFoundError(qualified_name)
        module = module_from_spec(spec)
        try:
            sys.modules[qualified_name] = module
            spec.loader.exec_module(module)
        except FileNotFoundError as exc:
            sys.modules.pop(qualified_name, None)
            msg = (
                f"Optional module '{qualified_name}' could not be loaded. "
                f"Ensure '{module_path}' is included in your environment."
            )
            raise ModuleNotFoundError(msg) from exc
        return module

    def _load_from_resource(resource: Any) -> ModuleType | None:
        with resources.as_file(resource) as module_path:
            try:
                return _load_from_path(Path(module_path))
            except ModuleNotFoundError:
                return None

    try:
        package_files = resources.files(package_name)
    except ModuleNotFoundError:
        package_files = None

    if package_files is not None:
        for suffix in suffixes:
            resource_name = f"{module_name}{suffix}"
            module_resource = package_files.joinpath(resource_name)
            if module_resource.is_file():
                module = _load_from_resource(module_resource)
                if module is not None:
                    return module

    candidate_roots = _candidate_roots()
    for root in candidate_roots:
        for suffix in suffixes:
            candidate = root / f"{module_name}{suffix}"
            if candidate.is_file():
                try:
                    return _load_from_path(candidate)
                except ModuleNotFoundError:
                    continue

        package_init_candidates = [
            root / module_name / "__init__.py",
            root / module_name / "__init__.pyc",
        ]
        for package_candidate in package_init_candidates:
            if package_candidate.is_file():
                try:
                    return _load_from_path(package_candidate)
                except ModuleNotFoundError:
                    continue

        pycache = root / "__pycache__"
        if pycache.is_dir():
            pattern = f"{module_name}.*.pyc"
            for compiled in pycache.glob(pattern):
                try:
                    return _load_from_path(compiled)
                except ModuleNotFoundError:
                    continue

    resource_list = ", ".join(f"{module_name}{suffix}" for suffix in suffixes)
    msg = (
        f"Optional module '{qualified_name}' is not available. "
        f"Expected one of ({resource_list}) to be present in this installation. "
        f"Checked directories: {[str(root) for root in candidate_roots]}"
    )
    raise ModuleNotFoundError(msg)


def _export_from_module(module: ModuleType, names: tuple[str, ...]) -> None:
    """Populate the module globals with the selected attributes."""

    for name in names:
        globals()[name] = getattr(module, name)


def _import_optional(submodule: str) -> ModuleType:
    """Import ``submodule`` from this package with a local fallback."""

    qualified_name = f"{__name__}.{submodule}"
    try:
        return import_module(qualified_name)
    except ModuleNotFoundError as exc:
        if exc.name != qualified_name:
            raise
    return _load_local_module(submodule)


def _collect_exports(
    *groups: tuple[str, ...], extras: tuple[str, ...] = (),
) -> list[str]:
    """Collect export names from *groups* while preserving order and uniqueness."""

    seen: set[str] = set()
    ordered: list[str] = []

    for name in itertools.chain.from_iterable(groups + (extras,)):
        if name not in seen:
            seen.add(name)
            ordered.append(name)

    return ordered


catalog_module_name = f"{__name__}.catalog"
catalog_spec = find_spec(catalog_module_name)

if catalog_spec is not None:
    catalog_module = import_module(catalog_module_name)
else:  # pragma: no cover - environment-specific fallback
    catalog_module = _load_local_module("catalog")

_export_from_module(catalog_module, _CATALOG_EXPORTS)


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

cli_module = _import_optional("cli")
_export_from_module(cli_module, _CLI_EXPORTS)

pubchem_module = _import_optional("pubchem")
_export_from_module(pubchem_module, _PUBCHEM_EXPORTS)
from library.clients import pubchem as pc
from library.common.csv_utils import (
    write_csv_chunks_deterministic as write_csv_deterministic,
)
from library.common.log import logger
from library.integration import pubchem_library as pl
from library.integration.chembl_client import ChemblClient
from library.metadata import file_sha256, write_meta_yaml
from library.table_quality import analyze_table_quality
from library.validation import validate_testitems

_LAZY_EXPORTS: dict[str, str] = {
    "ReadInputIdsResult": "library.testitem_pipeline.cli:ReadInputIdsResult",
    "TestitemFetchError": "library.testitem_pipeline.cli:TestitemFetchError",
    "TestitemPipelineOptions": "library.testitem_pipeline.cli:TestitemPipelineOptions",
    "TestitemPipelineStageError": "library.testitem_pipeline.cli:TestitemPipelineStageError",
    "_FETCH_ERROR_SAMPLE_SIZE": "library.testitem_pipeline.cli:_FETCH_ERROR_SAMPLE_SIZE",
    "_log_missing_identifier_summary": "library.testitem_pipeline.cli:_log_missing_identifier_summary",
    "_prepare_pubchem_api_cfg": "library.testitem_pipeline.cli:_prepare_pubchem_api_cfg",
    "apply_testitem_enrichment": "library.testitem_pipeline.cli:apply_testitem_enrichment",
    "fetch_testitems": "library.testitem_pipeline.cli:fetch_testitems",
    "finalize_output": "library.testitem_pipeline.cli:finalize_output",
    "integrate_missing_identifiers": "library.testitem_pipeline.cli:integrate_missing_identifiers",
    "read_input_ids": "library.testitem_pipeline.cli:read_input_ids",
    "run_testitem_pipeline": "library.testitem_pipeline.cli:run_testitem_pipeline",
    "PUBCHEM_CID_CACHE_ENCODING": "library.testitem_pipeline.pubchem:PUBCHEM_CID_CACHE_ENCODING",
    "PUBCHEM_COLUMNS": "library.testitem_pipeline.pubchem:PUBCHEM_COLUMNS",
    "_CID_CACHE_MISSING": "library.testitem_pipeline.pubchem:_CID_CACHE_MISSING",
    "_PUBCHEM_CACHE_SCHEMA_VERSION": "library.testitem_pipeline.pubchem:_PUBCHEM_CACHE_SCHEMA_VERSION",
    "_PUBCHEM_SESSION_LOCK": "library.testitem_pipeline.pubchem:_PUBCHEM_SESSION_LOCK",
    "_PUBCHEM_SESSION_SIGNATURE": "library.testitem_pipeline.pubchem:_PUBCHEM_SESSION_SIGNATURE",
    "_load_pubchem_cid_cache": "library.testitem_pipeline.pubchem:_load_pubchem_cid_cache",
    "_merge_pubchem_properties": "library.testitem_pipeline.pubchem:_merge_pubchem_properties",
    "_prepare_pubchem_caches": "library.testitem_pipeline.pubchem:_prepare_pubchem_caches",
    "_prefetch_parents": "library.testitem_pipeline.pubchem:_prefetch_parents",
    "_pubchem_session_signature": "library.testitem_pipeline.pubchem:_pubchem_session_signature",
    "_resolve_pubchem_cids": "library.testitem_pipeline.pubchem:_resolve_pubchem_cids",
    "_write_pubchem_cid_cache": "library.testitem_pipeline.pubchem:_write_pubchem_cid_cache",
    "add_pubchem_data": "library.testitem_pipeline.pubchem:add_pubchem_data",
    "augment_pubchem": "library.testitem_pipeline.pubchem:augment_pubchem",
    "resolve_pubchem_cid": "library.testitem_pipeline.pubchem:resolve_pubchem_cid",
    "pl": "library.integration.pubchem_library",
}


def __getattr__(name: str) -> Any:
    """Динамически импортировать CLI и PubChem-хелперы при первом обращении."""

    if name in _LAZY_EXPORTS:
        module_name, _, attribute = _LAZY_EXPORTS[name].partition(":")
        module = import_module(module_name)
        value = getattr(module, attribute) if attribute else module
        globals()[name] = value
        return value
    raise AttributeError(f"module '{__name__}' не имеет атрибута '{name}'")


__all__ = [
    "ChemblClient",
    "analyze_table_quality",
    "file_sha256",
    "json",
    "logger",
    "pc",
    "pl",
    "requests",
    "testitem_enrichment",
    "validate_testitems",
    "write_csv_deterministic",
    "write_parent_catalog_cache",
    "write_meta_yaml",
]

_COLLECTED_EXTRAS = globals().get("_EXTRA_EXPORTS", ())

__all__ = _collect_exports(
    _CATALOG_EXPORTS,
    _CLI_EXPORTS,
    _PUBCHEM_EXPORTS,
    extras=_COLLECTED_EXTRAS,
)
