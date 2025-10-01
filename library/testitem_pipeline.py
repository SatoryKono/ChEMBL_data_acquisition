"""Reusable components for the test item data acquisition pipeline."""

from __future__ import annotations

import json
import sys
from collections import ChainMap
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Mapping, MutableMapping, NamedTuple, Sequence

import pandas as pd
import requests
from pandera.errors import SchemaErrors

from library import chembl_library as cl
from library import io, molecule_catalog, testitem_enrichment
from library import pubchem_library as pl
from library.chembl_client import ChemblClient
from library.config import (
    ApiCfg,
    Config,
    IoCfg,
    MoleculeCatalogCfg,
    PubChemCfg,
    RetryCfg,
    _serialize_paths,
)
from library.csv_utils import write_csv_deterministic
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.molecule_catalog import (
    load_parent_catalog,
    query_parent_catalog,
    update_parent_catalog_cache,
    write_parent_catalog_cache,
)
from library.pipeline_metadata import add_pipeline_metadata
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from library.validation import validate_testitems
from library.utils.atomic import open_atomic
from schemas import TestitemsSchema, normalize_testitems


UTC = timezone.utc  # noqa: UP017
ENCODING_UTF8 = "utf-8"
PUBCHEM_CID_CACHE_ENCODING = ENCODING_UTF8
_TYPO_PARENT_COLUMN = "parant_molecule_id"
_MOLECULE_HIERARCHY_COLUMNS = (
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
)

PUBCHEM_COLUMNS = [
    "pubchem_cid",
    "pubchem_iupac_name",
    "pubchem_molecular_formula",
    "pubchem_isomeric_smiles",
    "pubchem_canonical_smiles",
    "pubchem_inchi",
    "pubchem_inchikey",
]

_CID_CACHE_MISSING = object()
_PUBCHEM_CACHE_SCHEMA_VERSION = 1

_DEFAULT_CATALOG_CFG = MoleculeCatalogCfg()


def _pubchem_session_signature(api_cfg: ApiCfg, retry_cfg: RetryCfg) -> str:
    """Return a stable signature for the PubChem session configuration."""

    payload = {
        "api": api_cfg.model_dump(mode="json"),
        "retry": retry_cfg.model_dump(mode="json"),
    }
    return json.dumps(payload, sort_keys=True)


_PUBCHEM_SESSION_SIGNATURE: str | None = None


@dataclass
class ParentEnrichmentPreparation:
    """Intermediate data required to attach parent identifiers."""

    df: pd.DataFrame
    lookup_data: "ParentLookupPreparedData"
    parent_catalog: dict[str, str] | None
    parent_catalog_source: str
    parent_stats: "ParentLookupStats"


@dataclass
class ParentEnrichmentResult:
    """Result returned after running the parent enrichment stage."""

    df: pd.DataFrame
    parent_stats: "ParentLookupStats"


@dataclass(frozen=True)
class ParentLookupStats:
    """Summary information about parent molecule enrichment."""

    source: str
    missing: int
    unique: int
    attached: int
    uncovered: int


class ParentLookupPreparedData(NamedTuple):
    """Container for precomputed parent lookup data."""

    child_ids: pd.Series
    existing_parent_ids: pd.Series
    need_lookup: set[str]


PARENT_LOOKUP_SOURCE_CACHE = "cache"
PARENT_LOOKUP_SOURCE_LOOKUP = "lookup"
PARENT_LOOKUP_SOURCE_PARTIAL = "partial"
PARENT_LOOKUP_SOURCE_SYNC = "sync"
PARENT_LOOKUP_SOURCE_SKIPPED = "skipped"


def ensure_no_parant_column(df: pd.DataFrame) -> None:
    """Raise a :class:`ValueError` if the legacy typo column is present."""

    if _TYPO_PARENT_COLUMN in df.columns:
        raise ValueError(
            "unexpected column 'parant_molecule_id'; use 'parent_molecule_id' instead"
        )


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


@lru_cache(maxsize=None)
def _load_molecule_hierarchy_mapping(
    path: str,
    encoding: str,
    delimiter: str,
) -> dict[str, str | None]:
    """Return cached child → parent mapping with normalised identifiers."""

    csv_path = Path(path)
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Molecule hierarchy dictionary not found: {csv_path}"
        )

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

    missing_columns = [
        column for column in _MOLECULE_HIERARCHY_COLUMNS if column not in frame.columns
    ]
    if missing_columns:
        raise ValueError(
            "Molecule hierarchy dictionary missing required columns: "
            + ", ".join(missing_columns)
        )

    subset = frame.loc[:, list(_MOLECULE_HIERARCHY_COLUMNS)].copy()
    for column in _MOLECULE_HIERARCHY_COLUMNS:
        subset[column] = (
            subset[column]
            .fillna("")
            .astype("string")
            .str.strip()
            .str.upper()
        )

    subset = subset[subset["molecule_chembl_id"] != ""]
    subset = subset.drop_duplicates(
        subset=["molecule_chembl_id"],
        keep="first",
    )

    lookup: dict[str, str | None] = {}
    for molecule_id, parent_id in subset.itertuples(index=False, name=None):
        parent = parent_id or None
        if parent is not None and parent == molecule_id:
            parent = None
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
            "parent_molecule_chembl_id": value,
        }
        for key, value in cached.items()
    }


def load_molecule_hierarchy_lookup(
    path: Path | None,
    *,
    io_cfg: IoCfg,
    encoding: str | None = None,
    delimiter: str | None = None,
    catalog_cfg: MoleculeCatalogCfg | None = None,
) -> dict[str, str]:
    """Return child → parent mapping loaded from a local hierarchy file."""

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

    lookup = dict(raw_lookup)
    if not lookup:
        return {}

    attached_rows = sum(1 for value in lookup.values() if value is not None)

    logger.info(
        "molecule_hierarchy_lookup_loaded",
        path=str(resolved_path),
        rows=attached_rows,
    )
    return lookup


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

    if missing_ids:
        logger.warning(
            "parent_lookup_missing_parents",
            count=len(missing_ids),
            identifiers=missing_ids,
        )

    refreshed_parent = normalised_child.map(parent_map).astype("string")

    combined_parent = existing_parent.astype("string").copy()

    update_mask = combined_parent.isna() | combined_parent.eq("")
    combined_parent.loc[update_mask] = refreshed_parent.loc[update_mask]

    if getattr(catalog_cfg, "force_refresh_existing", False):
        existing_normalised = combined_parent.fillna("").astype("string")
        refreshed_normalised = refreshed_parent.fillna("").astype("string")
        mismatch_mask = existing_normalised != refreshed_normalised
        if mismatch_mask.any():
            combined_parent.loc[mismatch_mask] = refreshed_parent.loc[mismatch_mask]

    result[parent_column] = combined_parent.astype("string")

    missing = int(combined_parent.isna().sum())
    attached = int(combined_parent.notna().sum())

    final_source = source_resolved
    if catalog is not None and not missing_ids:
        if final_source in (None, PARENT_LOOKUP_SOURCE_CACHE, PARENT_LOOKUP_SOURCE_SKIPPED):
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
        missing=int(missing),
        unique=int(len(unique_children)),
        attached=int(attached),
        uncovered=int(uncovered_children),
    )

    logger.info(
        "parent_lookup_progress",
        source=stats.source,
        unique=stats.unique,
        attached=stats.attached,
        missing=stats.missing,
        uncovered=stats.uncovered,
    )

    return result, stats


def prepare_parent_enrichment(
    df: pd.DataFrame,
    *,
    catalog_cfg: MoleculeCatalogCfg,
    io_cfg: IoCfg,
    api_cfg: ApiCfg,
    timeout: float,
    client: ChemblClient,
    hierarchy_lookup_path: Path | None,
) -> tuple[int, ParentEnrichmentPreparation | None]:
    """Prepare DataFrame and catalog information prior to enrichment."""

    parent_stats = ParentLookupStats(
        source=PARENT_LOOKUP_SOURCE_SKIPPED,
        missing=0,
        unique=0,
        attached=0,
        uncovered=0,
    )

    parent_column = catalog_cfg.parent_field
    child_column = catalog_cfg.child_field

    if child_column in df.columns:
        normalised_ids = _normalise_chembl_ids(df[child_column])
    else:
        normalised_ids = pd.Series("", index=df.index, dtype="string")

    if parent_column in df.columns:
        existing_parent = _normalise_chembl_ids(df[parent_column])
    else:
        existing_parent = pd.Series("", index=df.index, dtype="string")

    resolved_lookup_path = hierarchy_lookup_path or catalog_cfg.hierarchy_lookup_path
    if resolved_lookup_path:
        try:
            hierarchy_lookup = load_molecule_hierarchy_lookup(
                resolved_lookup_path,
                io_cfg=io_cfg,
                encoding=catalog_cfg.hierarchy_lookup_encoding,
                delimiter=catalog_cfg.hierarchy_lookup_delimiter,
                catalog_cfg=catalog_cfg,
            )
        except ValueError as exc:
            logger.error(
                "molecule_hierarchy_lookup_invalid",
                error=str(exc),
                path=str(resolved_lookup_path),
            )
            return 1, None
    else:
        hierarchy_lookup = {}

    if hierarchy_lookup:
        hierarchy_values = {
            key: value for key, value in hierarchy_lookup.items() if value is not None
        }
    else:
        hierarchy_values = {}

    if hierarchy_values:
        hierarchy_series = normalised_ids.map(
            lambda value: hierarchy_values.get(value) if value else None
        )
        hierarchy_mask = hierarchy_series.notna()
        if hierarchy_mask.any():
            resolved = hierarchy_series[hierarchy_mask].astype("string")
            if parent_column in df.columns:
                df[parent_column] = df[parent_column].astype("string")
            else:
                df[parent_column] = pd.Series(pd.NA, index=df.index, dtype="string")
            df.loc[hierarchy_mask, parent_column] = resolved.astype(object)
            existing_parent.loc[hierarchy_mask] = resolved.fillna("").astype("string")

    if getattr(catalog_cfg, "force_refresh_existing", False):
        need_lookup_mask = normalised_ids != ""
    else:
        need_lookup_mask = (normalised_ids != "") & (existing_parent == "")
    initial_need_lookup = set(normalised_ids[need_lookup_mask])
    need_lookup = set(initial_need_lookup)

    cache_before = _cache_state(catalog_cfg.cache_path)
    cache_after = cache_before
    parent_catalog: dict[str, str] = {}
    parent_catalog_source = PARENT_LOOKUP_SOURCE_SKIPPED

    if need_lookup and cache_before[0]:
        try:
            parent_catalog = query_parent_catalog(
                need_lookup,
                catalog_cfg=catalog_cfg,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error(
                "parent_catalog_invalid",
                error=str(exc),
                path=str(catalog_cfg.cache_path),
            )
            return 1, None
        cache_after = _cache_state(catalog_cfg.cache_path)
        parent_catalog_source = (
            PARENT_LOOKUP_SOURCE_CACHE
            if cache_after == cache_before
            else PARENT_LOOKUP_SOURCE_SYNC
        )
        if parent_catalog:
            need_lookup -= set(parent_catalog)

    if need_lookup:
        try:
            fetched = molecule_catalog.fetch_parent_catalog_for(
                need_lookup,
                client=client,
                api_cfg=api_cfg,
                timeout=timeout,
                catalog_cfg=catalog_cfg,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error("parent_lookup_partial_fetch_failed", error=str(exc))
            return 1, None
        if fetched:
            parent_catalog.update(fetched)
            update_parent_catalog_cache(fetched, catalog_cfg)
            parent_catalog_source = PARENT_LOOKUP_SOURCE_PARTIAL
            need_lookup -= set(fetched)

    if need_lookup:
        try:
            fallback_catalog = load_parent_catalog(
                client=client,
                api_cfg=api_cfg,
                catalog_cfg=catalog_cfg,
                timeout=timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error("parent_catalog_invalid", error=str(exc))
            return 1, None
        if fallback_catalog:
            parent_catalog.update(fallback_catalog)
            parent_catalog_source = PARENT_LOOKUP_SOURCE_CACHE
            need_lookup -= set(fallback_catalog)

    lookup_resolved = initial_need_lookup - need_lookup
    parent_lookup_data = ParentLookupPreparedData(
        child_ids=normalised_ids,
        existing_parent_ids=existing_parent,
        need_lookup=set(need_lookup),
    )
    if (
        lookup_resolved
        and not need_lookup
        and parent_catalog_source
        not in (PARENT_LOOKUP_SOURCE_PARTIAL, PARENT_LOOKUP_SOURCE_SYNC)
    ):
        parent_catalog_source = PARENT_LOOKUP_SOURCE_LOOKUP

    try:
        ensure_no_parant_column(df)
    except ValueError as exc:
        logger.error(
            "invalid_column",
            column=_TYPO_PARENT_COLUMN,
            error=str(exc),
        )
        return 1, None

    return (
        0,
        ParentEnrichmentPreparation(
            df=df,
            lookup_data=parent_lookup_data,
            parent_catalog=parent_catalog or None,
            parent_catalog_source=parent_catalog_source,
            parent_stats=parent_stats,
        ),
    )


def run_parent_enrichment(
    prep: ParentEnrichmentPreparation,
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    catalog_cfg: MoleculeCatalogCfg,
    timeout: float,
) -> tuple[int, ParentEnrichmentResult | None]:
    """Attach parent molecule identifiers using the prepared context."""

    logger.info("parent_lookup_start")
    try:
        df, parent_stats = attach_parent_molecule_ids(
            prep.df,
            client=client,
            api_cfg=api_cfg,
            catalog_cfg=catalog_cfg,
            timeout=timeout,
            catalog=prep.parent_catalog,
            source=prep.parent_catalog_source,
            precomputed=prep.lookup_data,
        )
    except (requests.RequestException, ValueError) as exc:
        logger.error("parent_lookup_failed", error=str(exc))
        return 1, None

    logger.info(
        "parent_lookup_done",
        source=parent_stats.source,
        unique=parent_stats.unique,
        attached=parent_stats.attached,
        missing=parent_stats.missing,
        uncovered=parent_stats.uncovered,
    )

    return 0, ParentEnrichmentResult(df=df, parent_stats=parent_stats)


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
        primary = _select_primary_cid(
            value,
            chembl_id=key,
            identifier="cache_file",
            value=value,
        )
        if primary is not None:
            cache[key] = primary
    return cache


def _write_pubchem_cid_cache(
    path: Path | None, cache: Mapping[str, str | None]
) -> None:
    """Persist CID cache mapping to disk."""

    if path is None:
        return
    serialisable: dict[str, str] = {}
    for key, value in cache.items():
        if not key:
            continue
        if value is None:
            continue
        serialisable[key] = value
    try:
        with open_atomic(path, encoding=PUBCHEM_CID_CACHE_ENCODING) as handle:
            payload = {
                "metadata": {
                    "version": _PUBCHEM_CACHE_SCHEMA_VERSION,
                    "updated_at": datetime.now(UTC).isoformat(),
                },
                "values": serialisable,
            }
            json.dump(payload, handle, indent=2, sort_keys=True)
    except OSError as exc:  # pragma: no cover - I/O errors
        logger.warning("pubchem_cache_write_failed", path=str(path), error=str(exc))


def _select_primary_cid(
    candidates: str | None,
    *,
    chembl_id: str | None,
    identifier: str,
    value: str | None,
) -> str | None:
    """Return the primary CID from a pipe-separated list."""

    if not candidates:
        return None
    cid_list = [cid.strip() for cid in str(candidates).split("|") if cid.strip()]
    if not cid_list:
        return None
    primary = cid_list[0]
    if len(cid_list) > 1:
        logger.info(
            "pubchem_multiple_cid",
            chembl_id=chembl_id,
            identifier=identifier,
            value=value,
            cid=primary,
            alternatives=cid_list[1:],
        )
    return primary


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


def _pubchem_resolution_key(row: pd.Series) -> tuple[str | None, ...]:
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
    resolution_cache: MutableMapping[tuple[str | None, ...], pl.PubChemResolution]
    | None = None,
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

    parent_id = _normalise_identifier(
        row.get("parent_molecule_chembl_id"), uppercase=True
    )
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


def add_pubchem_data(
    df: pd.DataFrame,
    cfg: PubChemCfg,
    *,
    client: ChemblClient | None = None,
    api_cfg: ApiCfg | None = None,
    timeout: float | None = None,
    cid_cache: MutableMapping[str, str | None] | None = None,
    resolution_cache: MutableMapping[tuple[str | None, ...], pl.PubChemResolution]
    | None = None,
    parent_record_cache: MutableMapping[str, pd.Series | None] | None = None,
    testitem_fields: Sequence[str] | None = None,
    request_limit: int = 0,
) -> pd.DataFrame:
    """Return ``df`` augmented with PubChem annotations."""

    if df.empty:
        return df

    if not getattr(cfg, "enable", True):
        logger.info("pubchem_disabled")
        return df

    result = df.reset_index(drop=True).copy()
    prefer_local = getattr(cfg, "prefer_local_smiles", False)
    existing_cols = [col for col in PUBCHEM_COLUMNS if col in result.columns]
    if prefer_local and existing_cols:
        normalised = pd.DataFrame(
            {
                col: result[col].astype("string").replace("", pd.NA)
                for col in existing_cols
            }
        )
        complete_mask = normalised.notna().all(axis=1)
    else:
        complete_mask = pd.Series(False, index=result.index)

    def _mask_for(series: pd.Series | None, keyword: str) -> pd.Series:
        if series is None:
            return pd.Series(False, index=result.index)
        normalised = series.astype("string").fillna("").str.strip().str.casefold()
        return normalised.str.contains(keyword.casefold(), regex=False, na=False)

    molecule_type_series = result.get("molecule_type")
    structure_type_series = result.get("structure_type")
    polymer_mask = _mask_for(molecule_type_series, "polymer") | _mask_for(
        structure_type_series, "polymer"
    )
    mixture_mask = _mask_for(molecule_type_series, "mixture") | _mask_for(
        structure_type_series, "mixture"
    )

    polymer_indexes = set(polymer_mask[polymer_mask].index.tolist())
    mixture_indexes = set(mixture_mask[mixture_mask].index.tolist())

    skip_indexes: set[int] = set()
    allow_polymer = getattr(cfg, "allow_polymer", False)
    if not allow_polymer:
        skip_indexes = polymer_indexes | mixture_indexes
        if skip_indexes:
            logger.warning(
                "pubchem_skip_polymers",
                count=len(skip_indexes),
                polymer_count=len(skip_indexes & polymer_indexes),
                mixture_count=len(skip_indexes & mixture_indexes),
                indexes=[int(index) for index in sorted(skip_indexes)],
            )

    cache_path = getattr(cfg, "cid_cache_path", None)
    cache_ttl_hours = getattr(cfg, "cache_ttl_hours", None)
    if cid_cache is None:
        cid_cache = _load_pubchem_cid_cache(cache_path, ttl_hours=cache_ttl_hours)
    cache_dirty = False
    if resolution_cache is None:
        resolution_cache = {}

    if parent_record_cache is None:
        parent_record_cache = {}

    local_records: pd.DataFrame | None = None
    if "molecule_chembl_id" in result.columns:
        prepared_local_records = (
            result.assign(
                __local_molecule=lambda frame: frame[
                    "molecule_chembl_id"
                ]
                .astype("string")
                .str.strip()
                .str.upper()
            )
            .dropna(subset=["__local_molecule"])
        )
        if not prepared_local_records.empty:
            local_records = (
                prepared_local_records.drop_duplicates("__local_molecule")
                .set_index("__local_molecule")
                .rename_axis("molecule_chembl_id")
                .reindex(columns=result.columns)
            )
            for chembl_norm in local_records.index:
                try:
                    record = local_records.loc[chembl_norm]
                except KeyError:  # pragma: no cover - defensive
                    continue
                parent_record_cache[chembl_norm] = record.copy()

    pending_parent_ids: list[str] = []
    seen_parent_ids: set[str] = set()
    if "parent_molecule_chembl_id" in result.columns:
        for parent_value in result["parent_molecule_chembl_id"]:
            parent_norm = _normalise_identifier(parent_value, uppercase=True)
            if not parent_norm:
                continue
            if parent_norm in parent_record_cache or (
                local_records is not None and parent_norm in local_records.index
            ):
                continue
            if parent_norm in seen_parent_ids:
                continue
            seen_parent_ids.add(parent_norm)
            pending_parent_ids.append(parent_norm)

    if pending_parent_ids and client is not None and api_cfg is not None:
        logger.info(
            "pubchem_parent_prefetch",
            count=len(pending_parent_ids),
            batch_size=getattr(cfg, "batch_size", 0),
        )
        try:
            fetched = cl.get_testitem(
                pending_parent_ids,
                cfg=api_cfg,
                client=client,
                chunk_size=getattr(cfg, "batch_size", 0),
                timeout=timeout,
                fields=testitem_fields,
                page_limit=request_limit,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.warning("pubchem_parent_prefetch_failed", error=str(exc))
        else:
            if not fetched.empty:
                fetched = fetched.astype("string")
                fetched["molecule_chembl_id"] = fetched["molecule_chembl_id"].str.upper()
                fetched = fetched.drop_duplicates("molecule_chembl_id")
                fetched = fetched.set_index("molecule_chembl_id")
                for parent_norm, row in fetched.iterrows():
                    parent_record_cache[parent_norm] = row

    def load_parent_record(parent_norm: str) -> pd.Series | None:
        if parent_norm in parent_record_cache:
            return parent_record_cache[parent_norm]
        if local_records is not None and parent_norm in local_records.index:
            try:
                return local_records.loc[parent_norm]
            except KeyError:
                pass
        return parent_record_cache.get(parent_norm)

    if "molecule_chembl_id" in result.columns and "pubchem_cid" in result.columns:
        for chembl_raw, cid_raw in zip(
            result["molecule_chembl_id"], result["pubchem_cid"], strict=False
        ):
            chembl_id = _normalise_identifier(chembl_raw, uppercase=True)
            if not chembl_id:
                continue
            cid_value = _normalise_identifier(cid_raw)
            if not cid_value:
                continue
            before = cid_cache.get(chembl_id, _CID_CACHE_MISSING)
            cid_cache[chembl_id] = cid_value
            if before is _CID_CACHE_MISSING or before != cid_value:
                cache_dirty = True

    if "pubchem_cid" in result.columns:
        cid_series = result["pubchem_cid"].astype("string").copy()
    else:
        cid_series = pd.Series(pd.NA, index=result.index, dtype="string")

    skip_mask = pd.Series(False, index=result.index)
    if skip_indexes:
        skip_mask.loc[list(skip_indexes)] = True

    prefer_local_mask = (
        complete_mask.astype(bool)
        if prefer_local
        else pd.Series(False, index=result.index)
    )

    if "molecule_chembl_id" in result.columns:
        chembl_norm = result["molecule_chembl_id"].map(
            lambda value: _normalise_identifier(value, uppercase=True)
        )
    else:
        chembl_norm = pd.Series(
            [None] * len(result), index=result.index, dtype="object"
        )

    def _is_cached(chembl_id: str | None) -> bool:
        return bool(chembl_id and cid_cache.get(chembl_id))

    cached_mask = chembl_norm.map(_is_cached)
    needs_lookup_mask = (
        chembl_norm.notna()
        & ~skip_mask
        & ~prefer_local_mask
        & ~cached_mask
    )

    total = int(needs_lookup_mask.sum())
    if total:
        logger.info("pubchem_start", total=total)
    else:
        logger.info("pubchem_no_smiles")

    cid_series = cid_series.astype("string")
    lookup_cids: set[str] = set()
    for idx, chembl_id in chembl_norm[cached_mask].items():
        cached_value = cid_cache.get(chembl_id)
        if cached_value:
            cid_series.loc[idx] = cached_value
            lookup_cids.add(cached_value)

    for progress, row in enumerate(result.loc[needs_lookup_mask].itertuples(), start=1):
        logger.info("pubchem_progress", current=progress, total=total)
        idx = row.Index
        chembl_id = chembl_norm.loc[idx]
        before_present = bool(chembl_id and chembl_id in cid_cache)
        before_value = cid_cache[chembl_id] if before_present else _CID_CACHE_MISSING
        cid = resolve_pubchem_cid(
            result.loc[idx],
            cid_cache,
            cfg,
            parent_loader=load_parent_record,
            resolution_cache=resolution_cache,
        )
        if cid:
            cid_series.loc[idx] = cid
            lookup_cids.add(cid)
        elif getattr(cfg, "write_not_found_literal", False):
            cid_series.loc[idx] = "Not Found"
        else:
            cid_series.loc[idx] = pd.NA
        if chembl_id:
            after_present = chembl_id in cid_cache
            after_value = cid_cache[chembl_id] if after_present else _CID_CACHE_MISSING
            if before_present != after_present or after_value != before_value:
                cache_dirty = True

    if cache_dirty:
        _write_pubchem_cid_cache(cache_path, cid_cache)

    def _value_or_na(value: str | None) -> object:
        return value if value not in (None, "") else pd.NA

    properties: dict[str, pl.Properties] = {}
    for cid in sorted(lookup_cids):
        properties[cid] = pl.get_properties(cid, cfg)

    properties_records: dict[str, dict[str, object]] = {}
    for cid, props in properties.items():
        properties_records[cid] = {
            "pubchem_iupac_name": _value_or_na(props.IUPACName),
            "pubchem_molecular_formula": _value_or_na(props.MolecularFormula),
            "pubchem_isomeric_smiles": _value_or_na(props.iSMILES),
            "pubchem_canonical_smiles": _value_or_na(props.cSMILES),
            "pubchem_inchi": _value_or_na(props.InChI),
            "pubchem_inchikey": _value_or_na(props.InChIKey),
        }

    properties_df = pd.DataFrame.from_dict(properties_records, orient="index")
    pubchem_df = cid_series.to_frame("pubchem_cid").join(properties_df, on="pubchem_cid")
    pubchem_df = pubchem_df.reindex(result.index)

    preserve_mask = skip_mask | prefer_local_mask
    if preserve_mask.any():
        existing_columns = [
            column
            for column in PUBCHEM_COLUMNS
            if column in result.columns and column in pubchem_df.columns
        ]
        if existing_columns:
            original_existing = result[existing_columns].astype("string")
            intersecting_columns = [
                column for column in existing_columns if column in pubchem_df.columns
            ]
            if intersecting_columns:
                pubchem_df[intersecting_columns] = (
                    pubchem_df[intersecting_columns]
                    .astype("string")
                    .mask(preserve_mask, original_existing[intersecting_columns])
                )
            missing_columns = [
                column for column in existing_columns if column not in pubchem_df.columns
            ]
            if missing_columns:
                pubchem_df[missing_columns] = original_existing[missing_columns]

    pubchem_df = pubchem_df.convert_dtypes()

    prefer_values = getattr(cfg, "prefer_local_values", False)

    for col in PUBCHEM_COLUMNS:
        if col not in pubchem_df.columns:
            pubchem_df[col] = pd.Series(pd.NA, index=result.index, dtype="string")
        new_series = pubchem_df[col].astype("string")
        if col in result.columns:
            existing = result[col].astype("string")
            result[col] = existing
            if prefer_local:
                missing_mask = existing.isna() | existing.eq("")
                result.loc[missing_mask, col] = new_series[missing_mask]
            elif prefer_values:
                incoming = new_series.where(
                    ~(new_series.isna() | new_series.eq("")), pd.NA
                )
                combined = incoming.combine_first(existing)
                result[col] = combined.astype("string")
            else:
                result[col] = new_series
        else:
            result[col] = new_series

    return result


def augment_pubchem(
    df: pd.DataFrame,
    *,
    pubchem_cfg: PubChemCfg,
    api_cfg: ApiCfg,
    retry_cfg: RetryCfg,
    timeout: float,
    client: ChemblClient,
    fields: Sequence[str] | None,
    request_limit: int,
) -> pd.DataFrame:
    """Augment ``df`` with PubChem information if enabled.

    The shared PubChem client session is initialised with ``api_cfg`` and
    ``retry_cfg`` before enrichment. Callers must therefore provide both
    configurations even when invoking this function directly.
    """

    pubchem_cid_cache: dict[str, str | None] | None = None
    pubchem_resolution_cache: (
        dict[tuple[str | None, ...], pl.PubChemResolution] | None
    ) = None
    pubchem_parent_record_cache: dict[str, pd.Series | None] | None = None
    if getattr(pubchem_cfg, "enable", True):
        global _PUBCHEM_SESSION_SIGNATURE
        session_signature = _pubchem_session_signature(api_cfg, retry_cfg)
        if session_signature != _PUBCHEM_SESSION_SIGNATURE:
            pl.init_session(api_cfg, retry_cfg)
            _PUBCHEM_SESSION_SIGNATURE = session_signature
        pubchem_cid_cache = _load_pubchem_cid_cache(
            getattr(pubchem_cfg, "cid_cache_path", None),
            ttl_hours=getattr(pubchem_cfg, "cache_ttl_hours", None),
        )
        pubchem_resolution_cache = {}
        pubchem_parent_record_cache = {}

    logger.info("pubchem_augment_start")
    result = add_pubchem_data(
        df,
        pubchem_cfg,
        client=client,
        api_cfg=api_cfg,
        timeout=timeout,
        cid_cache=pubchem_cid_cache,
        resolution_cache=pubchem_resolution_cache,
        parent_record_cache=pubchem_parent_record_cache,
        testitem_fields=fields,
        request_limit=request_limit,
    )
    logger.info("pubchem_augment_done")
    return result


def apply_testitem_enrichment(
    df: pd.DataFrame,
    *,
    enrichment_cfg,
    io_cfg: IoCfg,
) -> tuple[int, pd.DataFrame | None]:
    """Apply optional test item enrichment if enabled."""

    if not enrichment_cfg.enable:
        return 0, df

    logger.info("testitem_enrichment_start")
    try:
        enriched = testitem_enrichment.enrich(
            df,
            cfg=enrichment_cfg,
            io_cfg=io_cfg,
        )
    except ValueError as exc:
        logger.error("testitem_enrichment_failed", error=str(exc))
        return 1, None
    logger.info("testitem_enrichment_done")
    return 0, enriched


def finalize_output(
    df: pd.DataFrame,
    *,
    cfg: Config,
    output: Path,
    parent_stats: ParentLookupStats,
    input_csv: Path,
    rows_total: int,
    missing_ids: Sequence[str] | None = None,
) -> int:
    """Normalise, validate, and persist the final dataset."""

    df = normalize_testitems(df)
    if "pubchem_cid" in df.columns:
        df["pubchem_cid"] = df["pubchem_cid"].astype(object)
    df = add_pipeline_metadata(df)

    schema_cols = list(TestitemsSchema.columns)
    head = [c for c in schema_cols if c in df.columns]
    tail = sorted(c for c in df.columns if c not in schema_cols)
    col_order = head + tail

    exit_code = 0
    required_cols = {name for name, col in TestitemsSchema.columns.items() if col.required}
    optional_cols = set(TestitemsSchema.columns) - required_cols
    missing_required = required_cols - set(df.columns)
    missing_optional = optional_cols - set(df.columns)
    if not missing_required:
        if missing_optional:
            logger.warning(
                "optional_columns_missing",
                columns=sorted(missing_optional),
            )
        try:
            validation_result = validate_testitems(df, return_result=True)
        except SchemaErrors as exc:
            failure_path = Path(output).with_name(
                f"{Path(output).stem}_failure_cases.csv"
            )
            errors = SidecarErrors()
            for row in exc.failure_cases.to_dict("records"):
                errors.add_error(row)
            errors.save(failure_path)
            logger.error(
                "validation_failed",
                failures=len(exc.failure_cases),
                path=str(failure_path),
            )
            df = getattr(exc, "validated_data", df)
            exit_code = 1
        else:
            df = validation_result.data
            if not validation_result.failure_cases.empty:
                failure_path = Path(output).with_name(
                    f"{Path(output).stem}_failure_cases.csv"
                )
                errors = SidecarErrors()
                for row in validation_result.failure_cases.to_dict("records"):
                    errors.add_error(row)
                errors.save(failure_path)
                logger.error(
                    "validation_failed",
                    failures=len(validation_result.failure_cases),
                    path=str(failure_path),
                )
                df = validation_result.data
                exit_code = 1
    else:
        logger.warning(
            "validation_skipped",
            missing_columns=sorted(missing_required),
        )

    rows_kept = len(df)
    rows_dropped = rows_total - rows_kept

    try:
        key_cols = ["molecule_chembl_id"]
        csv_chunksize = cfg.io.csv_chunksize
        csv_path = write_csv_deterministic(
            df,
            output,
            key_cols=key_cols or None,
            col_order=col_order,
            chunksize=csv_chunksize,
            sort_chunksize=csv_chunksize,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
        )
        logger.info("write_done", rows=rows_kept, path=str(csv_path))
    except OSError as exc:
        logger.error(
            "write_fail",
            error=str(exc),
            path=str(output),
        )
        return 1

    missing_ids = tuple(missing_ids or ())

    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_kept,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(csv_path),
        "parent_lookup_source": parent_stats.source,
        "parent_lookup_missing": parent_stats.missing,
    }
    if missing_ids:
        stats["missing_molecule_ids"] = list(missing_ids)
        stats["missing_molecule_ids_count"] = len(missing_ids)
    write_meta_yaml(
        csv_path=csv_path,
        command=" ".join(sys.argv),
        config_subset=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(input_csv)},
        stats=stats,
        schema="TestitemsSchema",
    )
    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error(
            "quality_report_failed",
            error=str(exc),
            path=str(output),
        )
        return 1

    return exit_code


__all__ = [
    "ParentEnrichmentPreparation",
    "ParentEnrichmentResult",
    "ParentLookupStats",
    "ParentLookupPreparedData",
    "PARENT_LOOKUP_SOURCE_CACHE",
    "PARENT_LOOKUP_SOURCE_LOOKUP",
    "PARENT_LOOKUP_SOURCE_PARTIAL",
    "PARENT_LOOKUP_SOURCE_SYNC",
    "PARENT_LOOKUP_SOURCE_SKIPPED",
    "LoadMoleculeHierarchyLookup",
    "load_molecule_hierarchy_lookup",
    "attach_parent_molecule_ids",
    "prepare_parent_enrichment",
    "run_parent_enrichment",
    "resolve_pubchem_cid",
    "add_pubchem_data",
    "augment_pubchem",
    "apply_testitem_enrichment",
    "finalize_output",
    "ensure_no_parant_column",
    "PUBCHEM_COLUMNS",
]
