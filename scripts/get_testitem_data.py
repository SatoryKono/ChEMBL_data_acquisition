"""Command line interface for retrieving ChEMBL test item data.

The module wraps :func:`library.testitem_pipeline.run_testitem_pipeline` while
exposing helpers that tests can import directly. Entry points return numeric
exit codes rather than terminating the interpreter to simplify orchestration.
The :func:`ensure_no_parant_column` helper guards against legacy CSV exports
that still include the misspelled parent identifier column.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import ChainMap
from dataclasses import dataclass
from functools import lru_cache
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import (
    Any,
    Callable,
    Iterable,
    Mapping,
    MutableMapping,
    NamedTuple,
    Sequence,
    Hashable,
    cast,
)

import pandas as pd
import requests


# ===== Parameters =====

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_NAME = "testitem.csv"
DEFAULT_OUTPUT_STEM = "testitems"


if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from library.utils.bootstrap import ensure_project_root


if __package__ in {None, ""}:
    ensure_project_root()

from library import cli  # noqa: F401 - re-exported for monkeypatching in tests
from library import io
from library.integration import molecule_catalog
from library.integration import pubchem_library as pl
from library.cli import LoggerConfig
from library.cli import build_parser as base_parser
from library.cli_utils import run_cli_command
from library.config import (
    ApiCfg,
    Config,
    IoCfg,
    MoleculeCatalogCfg,
    PubChemCfg,
)
from library.common.log import logger
from library.clients import pubchem as pc  # noqa: F401 - patched in tests
import library.testitem_pipeline as pipeline
from library.integration.chembl_client import ChemblClient
from library.testitem_pipeline import (
    PUBCHEM_CID_CACHE_ENCODING,
    PUBCHEM_COLUMNS,
    ReadInputIdsResult,
    TestitemPipelineOptions,
    _DEFAULT_CATALOG_CFG,
    _FETCH_ERROR_SAMPLE_SIZE,
    _MOLECULE_HIERARCHY_COLUMNS,
    _PUBCHEM_CACHE_SCHEMA_VERSION,
    _TYPO_PARENT_COLUMN,
    analyze_table_quality,
    ensure_no_parant_column,
    file_sha256,
    fetch_testitems,
    integrate_missing_identifiers,
    load_parent_catalog,
    query_parent_catalog,
    read_input_ids,
    run_testitem_pipeline,
    update_parent_catalog_cache,
    write_meta_yaml,
    write_parent_catalog_cache,
    _prepare_pubchem_api_cfg,
    _write_pubchem_cid_cache,
    PARENT_LOOKUP_SOURCE_CACHE,
    PARENT_LOOKUP_SOURCE_LOOKUP,
    PARENT_LOOKUP_SOURCE_PARTIAL,
    PARENT_LOOKUP_SOURCE_SKIPPED,
    PARENT_LOOKUP_SOURCE_SYNC,
)

_NO_PARENT_MARKERS = {"", "NULL", "NO PARENT"}

def _normalise_identifier(value: Any, *, uppercase: bool = False) -> str | None:
    """Return ``value`` normalised for PubChem lookup."""

    if value is None:
        return None
    if isinstance(value, str):
        normalised = value.strip()
    else:
        if pd.isna(value):
            return None
        normalised = str(value).strip()
    if not normalised:
        return None
    return normalised.upper() if uppercase else normalised


def _normalise_parent_identifier(value: object, *, child_id: str) -> str | None:
    """Return normalised parent identifier or ``None`` for missing markers."""

    if value is None or pd.isna(value):
        return None

    normalised_parent = str(value).strip().upper()
    if not normalised_parent:
        return None
    if normalised_parent in _NO_PARENT_MARKERS:
        return None
    if normalised_parent == child_id:
        return None
    return normalised_parent


@lru_cache(maxsize=None)
def _load_molecule_hierarchy_mapping(
    path: str,
    encoding: str,
    delimiter: str,
) -> dict[str, str | None]:
    """Return cached child → parent mapping with normalised identifiers."""

    csv_path = Path(path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Molecule hierarchy dictionary not found: {csv_path}")

    try:
        frame = pd.read_csv(
            csv_path,
            sep=delimiter,
            encoding=encoding,
            dtype="string",
        )
    except ValueError as exc:  # pragma: no cover - pandas raises on missing columns
        raise ValueError(
            "Unable to read molecule hierarchy dictionary; verify the CSV format."
        ) from exc

    child_column, parent_column = _MOLECULE_HIERARCHY_COLUMNS
    child_missing = child_column not in frame.columns
    parent_missing = parent_column not in frame.columns

    if child_missing:
        raise ValueError(
            "Molecule hierarchy dictionary missing required columns: "
            + child_column
        )

    if parent_missing:
        logger.warning(
            "molecule_hierarchy_missing_parent_column",
            column=parent_column,
            path=str(csv_path),
        )

    subset = frame.loc[:, [child_column]].copy()
    if parent_missing:
        subset[parent_column] = pd.Series(pd.NA, index=subset.index, dtype="string")
    else:
        subset[parent_column] = frame[parent_column].copy()

    for column in _MOLECULE_HIERARCHY_COLUMNS:
        subset[column] = subset[column].astype("string").str.strip().str.upper()

    subset = subset[subset["molecule_chembl_id"].notna()]
    subset = subset[subset["molecule_chembl_id"] != ""]
    subset = subset.drop_duplicates(
        subset=["molecule_chembl_id"],
        keep="first",
    )

    lookup: dict[str, str | None] = {}
    for molecule_id, parent_id in subset.itertuples(index=False, name=None):
        parent = _normalise_parent_identifier(parent_id, child_id=molecule_id)
        lookup[molecule_id] = parent

    return lookup


def LoadMoleculeHierarchyLookup(
    path: Path | str | None = None,
    *,
    encoding: str | None = None,
    delimiter: str | None = None,
    catalog_cfg: MoleculeCatalogCfg | None = None,
) -> dict[str, dict[str, str | None]]:
    """Return cached molecule hierarchy lookup keyed by ``molecule_chembl_id``."""

    cfg_source = catalog_cfg or _DEFAULT_CATALOG_CFG
    resolved_path_value = path or cfg_source.hierarchy_lookup_path
    if resolved_path_value is None:
        raise ValueError("hierarchy lookup path must be provided")
    resolved_path = Path(resolved_path_value)
    resolved_encoding = encoding or cfg_source.hierarchy_lookup_encoding
    resolved_delimiter = delimiter or cfg_source.hierarchy_lookup_delimiter

    cached = _load_molecule_hierarchy_mapping(
        str(resolved_path),
        resolved_encoding,
        resolved_delimiter,
    )
    return {
        key: {
            "molecule_chembl_id": key,
            "parent_molecule_chembl_id": _normalise_parent_identifier(
                value,
                child_id=key,
            ),
        }
        for key, value in cached.items()
    }


def _load_pubchem_cid_cache(
    path: Path | None, *, ttl_hours: float | None = None
) -> dict[str, str | None]:
    """Load CID cache mapping from disk."""

    if path is None:
        return {}
    try:
        with path.open("r", encoding=PUBCHEM_CID_CACHE_ENCODING) as handle:
            data = json.load(handle)
    except FileNotFoundError:
        return {}
    except (OSError, json.JSONDecodeError) as exc:  # pragma: no cover - I/O errors
        logger.warning("pubchem_cache_load_failed", path=str(path), error=str(exc))
        return {}
    metadata: Mapping[str, Any] | None = None
    values_source: Mapping[str, Any] | None = None
    if isinstance(data, Mapping) and "values" in data:
        values_candidate = data.get("values")
        if isinstance(values_candidate, Mapping):
            values_source = values_candidate
            metadata_candidate = data.get("metadata")
            if isinstance(metadata_candidate, Mapping):
                metadata = metadata_candidate
    elif isinstance(data, Mapping):
        values_source = data
    else:
        logger.warning("pubchem_cache_invalid", path=str(path))
        return {}

    if ttl_hours and ttl_hours > 0 and metadata is not None:
        updated_raw = metadata.get("updated_at")
        if isinstance(updated_raw, str):
            try:
                updated_at = datetime.fromisoformat(updated_raw)
            except ValueError:
                updated_at = None
            if updated_at is not None and updated_at.tzinfo is None:
                updated_at = updated_at.replace(tzinfo=UTC)
            if updated_at is not None:
                expires_at = updated_at + timedelta(hours=ttl_hours)
                if expires_at <= datetime.now(UTC):
                    logger.info("pubchem_cache_expired", path=str(path))
                    return {}

    cache: dict[str, str | None] = {}
    for raw_key, raw_value in (values_source or {}).items():
        key = _normalise_identifier(raw_key, uppercase=True)
        if not key:
            continue
        if raw_value is None:
            cache[key] = None
            continue
        value = _normalise_identifier(raw_value)
        if not value:
            continue
        primary = pl.select_primary_cid(
            value,
            identifier="cache_file",
            value=value,
            context_id=key,
        )
        if primary is not None:
            cache[key] = primary
    return cache


def _pubchem_identifiers(row: pd.Series) -> dict[str, str | None]:
    """Return mapping of identifier names to normalised values."""

    return {
        "pubchem_cid": _normalise_identifier(row.get("pubchem_cid")),
        "canonical_smiles": _normalise_identifier(row.get("canonical_smiles")),
        "pubchem_canonical_smiles": _normalise_identifier(
            row.get("pubchem_canonical_smiles")
        ),
        "isomeric_smiles": _normalise_identifier(row.get("isomeric_smiles")),
        "pubchem_isomeric_smiles": _normalise_identifier(
            row.get("pubchem_isomeric_smiles")
        ),
        "standard_inchi_key": _normalise_identifier(
            row.get("standard_inchi_key"), uppercase=True
        ),
        "pubchem_inchikey": _normalise_identifier(
            row.get("pubchem_inchikey"), uppercase=True
        ),
        "standard_inchi": _normalise_identifier(row.get("standard_inchi")),
        "pubchem_inchi": _normalise_identifier(row.get("pubchem_inchi")),
        "pref_name": _normalise_identifier(row.get("pref_name")),
    }


def _pubchem_resolution_key(row: pd.Series) -> Hashable:
    """Return a hashable key describing identifiers used for PubChem lookup."""

    identifiers = _pubchem_identifiers(row)
    return (
        identifiers.get("canonical_smiles"),
        identifiers.get("pubchem_canonical_smiles"),
        identifiers.get("isomeric_smiles"),
        identifiers.get("pubchem_isomeric_smiles"),
        identifiers.get("standard_inchi_key"),
        identifiers.get("pubchem_inchikey"),
        identifiers.get("standard_inchi"),
        identifiers.get("pubchem_inchi"),
        identifiers.get("pref_name"),
    )


def resolve_pubchem_cid(
    row: pd.Series,
    cache: MutableMapping[str, str | None],
    cfg: PubChemCfg,
    *,
    parent_loader: Callable[[str], pd.Series | None] | None = None,
    resolution_cache: MutableMapping[Hashable, pl.PubChemResolution] | None = None,
    visited: set[str] | None = None,
) -> str | None:
    """Resolve PubChem CID for a ChEMBL record."""

    chembl_id = _normalise_identifier(row.get("molecule_chembl_id"), uppercase=True)
    identifiers = _pubchem_identifiers(row)
    resolution_key = _pubchem_resolution_key(row)

    if visited is None:
        visited = set()

    if chembl_id:
        if chembl_id in visited:
            logger.info(
                "pubchem_parent_structure_missing",
                child=chembl_id,
                parent=chembl_id,
                reason="parent_cycle_detected",
            )
            if chembl_id not in cache:
                cache[chembl_id] = None
            return None
        visited.add(chembl_id)

    resolution = pl.resolve_pubchem_record(
        identifiers,
        cfg,
        cid_cache=cache,
        cache_key=chembl_id,
        resolution_cache=resolution_cache,
        resolution_key=resolution_key,
    )
    cid = resolution.cid
    if cid is not None:
        return cid

    if not cfg.use_parent_for_salts:
        if chembl_id and chembl_id not in cache:
            cache[chembl_id] = None
        return None

    parent_raw = None
    if isinstance(row, pd.Series) and "parent_molecule_chembl_id" in row.index:
        parent_raw = row.get("parent_molecule_chembl_id")
    parent_id = _normalise_identifier(parent_raw, uppercase=True)
    if not parent_id:
        if chembl_id and chembl_id not in cache:
            cache[chembl_id] = None
        return None

    if parent_loader is None:
        if chembl_id and chembl_id not in cache:
            cache[chembl_id] = None
        return None

    if parent_id in visited:
        logger.info(
            "pubchem_parent_structure_missing",
            child=chembl_id,
            parent=parent_id,
            reason="parent_cycle_detected",
        )
        if parent_id not in cache:
            cache[parent_id] = None
        if chembl_id and chembl_id not in cache:
            cache[chembl_id] = None
        return None

    parent_row = parent_loader(parent_id)
    if parent_row is None:
        logger.info(
            "pubchem_parent_structure_missing",
            child=chembl_id,
            parent=parent_id,
            reason="parent_unavailable",
        )
        if chembl_id and chembl_id not in cache:
            cache[chembl_id] = None
        return None

    parent_cid = resolve_pubchem_cid(
        parent_row,
        cache,
        cfg,
        parent_loader=parent_loader,
        resolution_cache=resolution_cache,
        visited=visited,
    )
    if parent_cid:
        cache[parent_id] = parent_cid
        if chembl_id:
            cache[chembl_id] = parent_cid
        return parent_cid

    logger.info(
        "pubchem_parent_structure_missing",
        child=chembl_id,
        parent=parent_id,
        reason="structure_unresolved",
    )
    if parent_id and parent_id not in cache:
        cache[parent_id] = None
    if chembl_id and chembl_id not in cache:
        cache[chembl_id] = None
    return None


@dataclass(frozen=True)
class ParentLookupStats:
    """Summary information about parent molecule enrichment."""

    source: str
    missing: int
    unique: int
    attached: int
    uncovered: int
    failed_ids: tuple[str, ...] = ()
    hierarchy_attached: int = 0
    fallback_attached: int = 0
    no_parent: int = 0


def _cache_state(path: Path) -> tuple[bool, float | None]:
    """Return cache file presence and modification time."""

    try:
        stat = path.stat()
    except FileNotFoundError:
        return False, None
    return True, stat.st_mtime


def _resolve_catalog_load_source(
    before: tuple[bool, float | None], after: tuple[bool, float | None]
) -> str:
    """Determine the lookup source after invoking ``load_parent_catalog``."""

    before_exists, before_mtime = before
    after_exists, after_mtime = after

    if after_exists and before_exists and after_mtime == before_mtime:
        return PARENT_LOOKUP_SOURCE_CACHE
    if after_exists:
        return PARENT_LOOKUP_SOURCE_SYNC
    if before_exists and not after_exists:
        return PARENT_LOOKUP_SOURCE_SYNC
    return PARENT_LOOKUP_SOURCE_SYNC


def _normalise_chembl_ids(series: pd.Series) -> pd.Series:
    """Return ``series`` normalised to upper-case ChEMBL identifiers."""

    normalised = series.astype("string").fillna("").str.strip().str.upper()
    return normalised


def load_molecule_hierarchy_lookup(
    path: Path | None,
    *,
    io_cfg: IoCfg,
    encoding: str | None = None,
    delimiter: str | None = None,
    catalog_cfg: MoleculeCatalogCfg | None = None,
) -> dict[str, str | None]:
    """Return child → parent mapping loaded from a local hierarchy file.

    The returned mapping mirrors
    :func:`library.testitem_pipeline.load_molecule_hierarchy_lookup` where
    values may be ``None`` when no parent is listed in the hierarchy.
    """

    cfg_source = catalog_cfg or _DEFAULT_CATALOG_CFG
    resolved_path_value = path or cfg_source.hierarchy_lookup_path
    if resolved_path_value is None:
        return {}

    resolved_path = Path(resolved_path_value)
    resolved_encoding = (
        encoding
        or getattr(cfg_source, "hierarchy_lookup_encoding", None)
        or io_cfg.csv_encoding
    )
    resolved_delimiter = (
        delimiter
        or getattr(cfg_source, "hierarchy_lookup_delimiter", None)
        or io_cfg.csv_sep
    )

    try:
        raw_lookup = _load_molecule_hierarchy_mapping(
            str(resolved_path),
            resolved_encoding,
            resolved_delimiter,
        )
    except FileNotFoundError:
        logger.info("molecule_hierarchy_lookup_missing", path=str(resolved_path))
        return {}
    except ValueError as exc:
        raise ValueError(f"invalid hierarchy lookup: {exc}") from exc

    lookup = {
        child_id: _normalise_parent_identifier(parent_id, child_id=child_id)
        for child_id, parent_id in raw_lookup.items()
    }
    if not lookup:
        return {}

    attached_rows = sum(1 for value in lookup.values() if value is not None)

    logger.info(
        "molecule_hierarchy_lookup_loaded",
        path=str(resolved_path),
        rows=attached_rows,
    )
    return lookup


class ParentLookupPreparedData(NamedTuple):
    """Container for precomputed parent lookup data."""

    child_ids: pd.Series
    existing_parent_ids: pd.Series
    need_lookup: set[str]


def attach_parent_molecule_ids(
    df: pd.DataFrame,
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    catalog_cfg: MoleculeCatalogCfg,
    timeout: float | None,
    catalog: Mapping[str, str] | None = None,
    source: str | None = None,
    precomputed: ParentLookupPreparedData | None = None,
) -> tuple[pd.DataFrame, ParentLookupStats]:
    """Attach parent molecule identifiers using the ChEMBL catalogue.

    The function operates on the original ``df`` content and can enrich frames
    that already contain parent identifiers as well as frames without the
    ``catalog_cfg.parent_field`` column.
    """

    if "parant_molecule_id" in df.columns:
        raise ValueError("Input frame contains unexpected column 'parant_molecule_id'.")

    result = df.copy()

    if result.empty:
        if catalog_cfg.parent_field not in result.columns:
            result[catalog_cfg.parent_field] = pd.Series(
                pd.NA, index=result.index, dtype="string"
            )
        stats = ParentLookupStats(
            source=PARENT_LOOKUP_SOURCE_SKIPPED,
            missing=0,
            unique=0,
            attached=0,
            uncovered=0,
        )
        return result, stats

    child_column = catalog_cfg.child_field
    parent_column = catalog_cfg.parent_field

    if child_column not in result.columns:
        logger.warning("parent_lookup_missing_child_column", column=child_column)
        result[parent_column] = pd.Series(pd.NA, index=result.index, dtype="string")
        stats = ParentLookupStats(
            source=PARENT_LOOKUP_SOURCE_SKIPPED,
            missing=len(result),
            unique=0,
            attached=0,
            uncovered=len(result),
        )
        return result, stats

    source_resolved = source
    if precomputed is not None:
        normalised_child = (
            precomputed.child_ids.reindex(result.index, fill_value="")
            .astype("string")
            .copy()
        )
        existing_parent = (
            precomputed.existing_parent_ids.reindex(result.index)
            .astype("string")
            .copy()
        )
        lookup_mask = normalised_child.isin(precomputed.need_lookup)
    else:
        normalised_child = _normalise_chembl_ids(result[child_column])
        if parent_column in result.columns:
            existing_parent = result[parent_column].astype("string").copy()
        else:
            existing_parent = pd.Series(pd.NA, index=result.index, dtype="string")
        lookup_mask = (normalised_child != "") & (
            existing_parent.isna() | existing_parent.eq("")
        )

    existing_parent_before = existing_parent.copy()

    unique_children = tuple(normalised_child[lookup_mask].unique().tolist())
    catalog_data: MutableMapping[str, str]
    used_partial_cache = False
    needs_full_sync = False
    partial_fetch_used = False
    full_sync_used = False
    uncovered_children = 0

    if catalog is not None:
        base_view = {key: catalog[key] for key in unique_children if key in catalog}
        catalog_data = ChainMap({}, base_view)
    else:
        catalog_data = {}
        if unique_children:
            queried = query_parent_catalog(unique_children, catalog_cfg)
            if queried:
                catalog_data.update(queried)
                used_partial_cache = True
                if source_resolved is None:
                    source_resolved = PARENT_LOOKUP_SOURCE_CACHE
            else:
                sqlite_exists = catalog_cfg.sqlite_path.is_file()
                used_partial_cache = sqlite_exists
                if sqlite_exists:
                    if source_resolved is None:
                        source_resolved = PARENT_LOOKUP_SOURCE_CACHE

    parent_map = {
        key: catalog_data[key] for key in unique_children if key in catalog_data
    }
    missing_ids = [key for key in unique_children if key not in parent_map]
    uncovered_children = len(missing_ids)

    fetched: dict[str, str] = {}
    if missing_ids and catalog is None:
        try:
            fetched = molecule_catalog.fetch_parent_catalog_for(
                missing_ids,
                client=client,
                api_cfg=api_cfg,
                timeout=timeout,
                catalog_cfg=catalog_cfg,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.warning("parent_lookup_partial_fetch_failed", error=str(exc))
            fetched = {}
        if fetched:
            partial_fetch_used = True
            catalog_data.update(fetched)
            parent_map.update(fetched)
            missing_ids = [key for key in unique_children if key not in parent_map]
            uncovered_children = len(missing_ids)
            if used_partial_cache:
                update_parent_catalog_cache(fetched, catalog_cfg)
            else:
                write_parent_catalog_cache(catalog_data, catalog_cfg)
                used_partial_cache = True

    needs_full_sync = catalog is None and uncovered_children > 0

    if missing_ids and catalog is None and needs_full_sync:
        cache_before_load = _cache_state(catalog_cfg.cache_path)
        catalog_data = load_parent_catalog(
            client=client,
            api_cfg=api_cfg,
            catalog_cfg=catalog_cfg,
            timeout=timeout,
        )
        cache_after_load = _cache_state(catalog_cfg.cache_path)
        full_sync_used = True
        source_resolved = _resolve_catalog_load_source(
            cache_before_load, cache_after_load
        )
        if partial_fetch_used:
            catalog_data.update(fetched)
        parent_map = {
            key: catalog_data.get(key, parent_map.get(key, ""))
            for key in unique_children
            if key in catalog_data or key in parent_map
        }
        missing_ids = [key for key in unique_children if key not in parent_map]
        uncovered_children = len(missing_ids)

    refreshed_parent = normalised_child.map(parent_map).astype("string")
    refreshed_normalised = refreshed_parent.fillna("").astype("string")

    combined_parent = existing_parent.copy()
    update_mask = combined_parent.isna() | combined_parent.eq("")
    combined_parent.loc[update_mask] = refreshed_parent.loc[update_mask]

    if getattr(catalog_cfg, "force_refresh_existing", False):
        existing_normalised = combined_parent.fillna("").astype("string")
        mismatch_mask = existing_normalised != refreshed_normalised
        if mismatch_mask.any():
            combined_parent.loc[mismatch_mask] = refreshed_parent.loc[mismatch_mask]
        existing_normalised = combined_parent.fillna("").astype("string")

    result[parent_column] = combined_parent.astype("string")

    final_normalised = combined_parent.fillna("").astype("string")
    previous_normalised = existing_parent_before.fillna("").astype("string")

    fallback_mask = (final_normalised != "") & (
        (previous_normalised == "") | (previous_normalised != final_normalised)
    )
    fallback_attached = int(fallback_mask.sum())
    no_parent_count = int((final_normalised == "").sum())
    attached = int((final_normalised != "").sum())
    missing = no_parent_count

    final_source = source_resolved
    if catalog is not None and not missing_ids:
        if final_source in (
            None,
            PARENT_LOOKUP_SOURCE_CACHE,
            PARENT_LOOKUP_SOURCE_SKIPPED,
        ):
            final_source = PARENT_LOOKUP_SOURCE_LOOKUP
    if full_sync_used:
        final_source = PARENT_LOOKUP_SOURCE_SYNC
    elif partial_fetch_used:
        final_source = PARENT_LOOKUP_SOURCE_PARTIAL
    elif final_source is None:
        final_source = (
            PARENT_LOOKUP_SOURCE_CACHE
            if used_partial_cache or catalog is not None
            else PARENT_LOOKUP_SOURCE_SYNC
        )

    stats = ParentLookupStats(
        source=final_source,
        missing=missing,
        unique=int(len(unique_children)),
        attached=int(attached),
        uncovered=int(uncovered_children),
        failed_ids=tuple(missing_ids),
        fallback_attached=int(fallback_attached),
        no_parent=int(no_parent_count),
    )

    logger.info(
        "parent_lookup_progress",
        source=stats.source,
        unique=stats.unique,
        attached=stats.attached,
        missing=stats.missing,
        uncovered=stats.uncovered,
        hierarchy_attached=stats.hierarchy_attached,
        fallback_attached=stats.fallback_attached,
        no_parent=stats.no_parent,
    )

    return result, stats


def add_pubchem_data(
    df: pd.DataFrame,
    cfg: PubChemCfg,
    *,
    client: ChemblClient | None = None,
    api_cfg: ApiCfg | None = None,
    timeout: float | None = None,
    cid_cache: MutableMapping[str, str | None] | None = None,
    resolution_cache: MutableMapping[Hashable, pl.PubChemResolution] | None = None,
    parent_record_cache: MutableMapping[str, pd.Series | None] | None = None,
    testitem_fields: Sequence[str] | None = None,
    request_limit: int = 1000,
) -> pd.DataFrame:
    """Augment ChEMBL records with PubChem information.

    Delegates to :func:`library.testitem_pipeline.add_pubchem_data` while
    relaxing the ``resolution_cache`` type to align with
    :func:`library.integration.pubchem_library.resolve_pubchem_record`.
    """

    return pipeline.add_pubchem_data(
        df,
        cfg,
        client=client,
        api_cfg=api_cfg,
        timeout=timeout,
        cid_cache=cid_cache,
        resolution_cache=cast(
            MutableMapping[tuple[str | None, ...], pl.PubChemResolution] | None,
            resolution_cache,
        ),
        parent_record_cache=parent_record_cache,
        testitem_fields=testitem_fields,
        request_limit=request_limit,
    )


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Invoke the reusable test item pipeline with CLI parameters.

    Parameters
    ----------
    cfg : Config
        Application configuration containing ChEMBL, PubChem and IO settings.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate that identifier loading,
        network requests or CSV export failed inside the test item pipeline.
    """

    output_csv = getattr(args, "output_csv", None)
    options = TestitemPipelineOptions(
        input_csv=Path(args.input_csv),
        output_csv=Path(output_csv) if output_csv else None,
    )
    return run_testitem_pipeline(cfg, options)


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the test item pipeline handling ``--skip-existing`` semantics."""

    output_path = Path(
        args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    )
    args.output_csv = output_path
    if args.skip_existing and output_path.exists() and not args.force:
        logger.info("pipeline_skip_existing", output=str(output_path))
        return 0
    logger.info(
        "testitem_pipeline_start",
        input=str(args.input_csv),
        output=str(output_path),
        limit=getattr(cfg.testitem, "limit", None),
        offset=getattr(args, "offset", getattr(cfg.testitem, "offset", None)),
        batch_size=getattr(cfg.testitem, "batch_size", None),
        timeout=getattr(cfg.testitem, "timeout", None),
    )
    exit_code = run_chembl(cfg, args)
    if exit_code == 0:
        logger.info("testitem_pipeline_done", output=str(output_path))
    else:
        logger.error(
            "testitem_pipeline_failed",
            output=str(output_path),
            exit_code=exit_code,
        )
    return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        Parser configured with the common CLI options and the associated logging
        configuration used by :func:`main`.
    """

    parser, log_cfg = base_parser(
        "ChEMBL and PubChem compound data utilities",
        column="molecule_chembl_id",
        chunk_size=1000,
        size_option="--batch-size",
        size_dest="batch_size",
    )
    parser.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help=(
            "Maximum number of identifiers to process; use 0 to skip processing"
        ),
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults.

    Parameters
    ----------
    argv : Sequence[str] | None, optional
        Command-line arguments to parse. When ``None`` the values from
        :data:`sys.argv` are used.

    Returns
    -------
    int
        ``0`` when the pipeline finishes successfully, non-zero otherwise.

    Raises
    ------
    SystemExit
        Raised when invalid command-line options are provided.
    """

    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)

    cli.prepare_io_paths(
        args,
        input_default=DEFAULT_INPUT_NAME,
        output_stem=DEFAULT_OUTPUT_STEM,
    )

    if args.limit == 0:
        logger.info("pipeline_skip_limit", limit=args.limit)
        return 0

    if args.limit is not None and args.limit < 0:
        parser.error("--limit must be zero or a positive integer")
    if args.offset < 0:
        parser.error("--offset must be zero or a positive integer")

    mapping = {
        "timeout": "testitem.timeout",
        "column": "testitem.column",
        "batch_size": "testitem.batch_size",
        "limit": "testitem.limit",
        "offset": "testitem.offset",
    }
    return run_cli_command(
        args=args,
        parser=parser,
        log_cfg=log_cfg,
        mapping=mapping,
        run=run,
        logger=logger,
    )


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
