"""Molecule catalog utilities for the test item pipeline."""

from __future__ import annotations

from collections import ChainMap
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, MutableMapping, NamedTuple, Sequence

import pandas as pd
import requests

from library.integration import molecule_catalog
from library.integration.chembl_client import ChemblClient
from library.config import ApiCfg, IoCfg, MoleculeCatalogCfg
from library.common.log import logger
from library.integration.molecule_catalog import (
    load_parent_catalog,
    query_parent_catalog,
    update_parent_catalog_cache,
    write_parent_catalog_cache,
)

_DEFAULT_CATALOG_CFG = MoleculeCatalogCfg()

_TYPO_PARENT_COLUMN = "parant_molecule_id"
_MOLECULE_HIERARCHY_COLUMNS = (
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
)
_NO_PARENT_MARKERS = {"", "NULL", "NO PARENT"}

PARENT_LOOKUP_SOURCE_CACHE = "cache"
PARENT_LOOKUP_SOURCE_LOOKUP = "lookup"
PARENT_LOOKUP_SOURCE_PARTIAL = "partial"
PARENT_LOOKUP_SOURCE_SYNC = "sync"
PARENT_LOOKUP_SOURCE_SKIPPED = "skipped"

_PARENT_SOURCE_PRIORITY: Mapping[str, int] = {
    PARENT_LOOKUP_SOURCE_SKIPPED: 0,
    PARENT_LOOKUP_SOURCE_CACHE: 1,
    PARENT_LOOKUP_SOURCE_LOOKUP: 2,
    PARENT_LOOKUP_SOURCE_PARTIAL: 3,
    PARENT_LOOKUP_SOURCE_SYNC: 4,
}


def _summarise_identifiers(
    values: Sequence[str], *, limit: int = 20
) -> tuple[list[str], bool]:
    """Return a sorted, truncated list of identifiers for structured logging."""

    identifiers = sorted({str(value) for value in values})
    truncated = len(identifiers) > limit
    if truncated:
        identifiers = identifiers[:limit]
    return identifiers, truncated


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
    failed_ids: tuple[str, ...] = ()
    hierarchy_attached: int = 0
    fallback_attached: int = 0
    no_parent: int = 0


class ParentLookupPreparedData(NamedTuple):
    """Container for precomputed parent lookup data."""

    child_ids: pd.Series
    existing_parent_ids: pd.Series
    need_lookup: set[str]


def _merge_parent_stats(base: ParentLookupStats, update: ParentLookupStats) -> ParentLookupStats:
    """Combine two :class:`ParentLookupStats` instances."""

    if base.source not in _PARENT_SOURCE_PRIORITY:
        base_priority = -1
    else:
        base_priority = _PARENT_SOURCE_PRIORITY[base.source]
    update_priority = _PARENT_SOURCE_PRIORITY.get(update.source, -1)
    if update_priority >= base_priority:
        resolved_source = update.source
    else:
        resolved_source = base.source
    combined_failed = tuple(
        dict.fromkeys((*base.failed_ids, *update.failed_ids))
    )
    return ParentLookupStats(
        source=resolved_source,
        missing=base.missing + update.missing,
        unique=base.unique + update.unique,
        attached=base.attached + update.attached,
        uncovered=base.uncovered + update.uncovered,
        failed_ids=combined_failed,
        hierarchy_attached=base.hierarchy_attached + update.hierarchy_attached,
        fallback_attached=base.fallback_attached + update.fallback_attached,
        no_parent=base.no_parent + update.no_parent,
    )


def ensure_no_parant_column(df: pd.DataFrame) -> None:
    """Raise a :class:`ValueError` if the legacy typo column is present."""

    if _TYPO_PARENT_COLUMN in df.columns:
        raise ValueError(
            "unexpected column 'parant_molecule_id'; use 'parent_molecule_id' instead"
        )


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


def load_molecule_hierarchy_lookup(
    path: Path | None,
    *,
    io_cfg: IoCfg,
    encoding: str | None = None,
    delimiter: str | None = None,
    catalog_cfg: MoleculeCatalogCfg | None = None,
) -> dict[str, str | None]:
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

    lookup = {
        child_id: _normalise_parent_identifier(parent_id, child_id=child_id)
        for child_id, parent_id in raw_lookup.items()
    }
    if not lookup:
        return {}

    attached_rows = sum(1 for value in lookup.values() if value is not None and value != "")

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
    """Attach parent molecule identifiers using the ChEMBL catalogue."""

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

    raw_unique_children = tuple(normalised_child[lookup_mask].unique().tolist())
    unique_children = tuple(
        value for value in raw_unique_children if not pd.isna(value)
    )
    catalog_data: MutableMapping[str, str]
    used_partial_cache = False
    needs_full_sync = False
    partial_fetch_used = False
    full_sync_used = False
    uncovered_children = 0
    parentless_filtered = molecule_catalog._filters_exclude_parentless(catalog_cfg)
    json_cache_exists = catalog_cfg.cache_path.is_file()
    sqlite_exists = catalog_cfg.sqlite_path.is_file()

    from library.pipelines import testitem as pipeline_module

    load_catalog_fn = getattr(pipeline_module, "load_parent_catalog", load_parent_catalog)
    query_catalog_fn = getattr(pipeline_module, "query_parent_catalog", query_parent_catalog)
    update_cache_fn = getattr(
        pipeline_module, "update_parent_catalog_cache", update_parent_catalog_cache
    )
    write_cache_fn = getattr(
        pipeline_module, "write_parent_catalog_cache", write_parent_catalog_cache
    )

    if catalog is not None:
        base_view = {key: catalog[key] for key in unique_children if key in catalog}
        catalog_data = ChainMap({}, base_view)
    else:
        catalog_data = {}
        if unique_children:
            queried = query_catalog_fn(unique_children, catalog_cfg)
            if queried:
                catalog_data.update(queried)
                used_partial_cache = True
                if source_resolved is None:
                    source_resolved = PARENT_LOOKUP_SOURCE_CACHE
            else:
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
                update_cache_fn(fetched, catalog_cfg)
            else:
                write_cache_fn(catalog_data, catalog_cfg)
                used_partial_cache = True

    needs_full_sync = catalog is None and uncovered_children > 0

    skip_full_sync = (
        missing_ids
        and catalog is None
        and needs_full_sync
        and parentless_filtered
        and not (json_cache_exists or sqlite_exists)
    )

    if skip_full_sync:
        identifiers, truncated = _summarise_identifiers(missing_ids)
        logger.info(
            "parent_lookup_full_sync_skipped_parentless",
            count=len(missing_ids),
            identifiers=identifiers,
            truncated=truncated,
        )
        source_resolved = PARENT_LOOKUP_SOURCE_SKIPPED
    elif missing_ids and catalog is None and needs_full_sync:
        cache_before_load = _cache_state(catalog_cfg.cache_path)
        cache_after_load = cache_before_load
        try:
            loaded_catalog = load_catalog_fn(
                client=client,
                api_cfg=api_cfg,
                catalog_cfg=catalog_cfg,
                timeout=timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.warning("parent_catalog_full_sync_failed", error=str(exc))
        else:
            catalog_data = loaded_catalog
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
        identifiers, truncated = _summarise_identifiers(missing_ids)
        log_missing = logger.warning
        if skip_full_sync and parentless_filtered:
            log_missing = logger.info
        log_missing(
            "parent_lookup_missing_parents",
            count=len(missing_ids),
            identifiers=identifiers,
            truncated=truncated,
        )

    refreshed_parent = normalised_child.map(parent_map).astype("string")
    refreshed_normalised = refreshed_parent.fillna("").astype("string")

    combined_parent = existing_parent.astype("string").copy()

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
        missing=int(missing),
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

    hierarchy_attached = 0

    dictionary_resolved_children: set[str] = set()
    if hierarchy_lookup:
        missing_sentinel = object()
        hierarchy_series = normalised_ids.map(
            lambda value: hierarchy_lookup.get(value, missing_sentinel)
            if value
            else missing_sentinel
        )
        hierarchy_series = pd.Series(hierarchy_series, index=df.index, dtype="object")
        hierarchy_mask = hierarchy_series.ne(missing_sentinel)
        if hierarchy_mask.any():
            resolved_values = hierarchy_series[hierarchy_mask]
            resolved_series = pd.Series(resolved_values, dtype="object")
            resolved_series = resolved_series.where(
                ~resolved_series.isna(),
                None,
            )
            resolved_as_string = resolved_series.astype("string")
            dictionary_resolved_children.update(
                normalised_ids.loc[hierarchy_mask & normalised_ids.ne("")].tolist()
            )
            if parent_column in df.columns:
                df[parent_column] = df[parent_column].astype("string")
            else:
                df[parent_column] = pd.Series(pd.NA, index=df.index, dtype="string")
            df.loc[hierarchy_mask, parent_column] = resolved_as_string
            existing_parent.loc[hierarchy_mask] = resolved_as_string
            hierarchy_attached = int(hierarchy_mask.sum())

            has_parent_mask = resolved_as_string.notna()
            if has_parent_mask.any():
                parent_updates = resolved_as_string[has_parent_mask]
                df.loc[parent_updates.index, parent_column] = parent_updates
                existing_parent.loc[parent_updates.index] = parent_updates

    if getattr(catalog_cfg, "force_refresh_existing", False):
        need_lookup_mask = normalised_ids != ""
    else:
        need_lookup_mask = ((normalised_ids != "") & (existing_parent == "")).fillna(False)
    initial_need_lookup = set(normalised_ids[need_lookup_mask])
    if dictionary_resolved_children:
        need_lookup = initial_need_lookup - dictionary_resolved_children
    else:
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

    parent_lookup_data = ParentLookupPreparedData(
        child_ids=normalised_ids,
        existing_parent_ids=existing_parent,
        need_lookup=need_lookup,
    )

    parent_stats = ParentLookupStats(
        source=parent_catalog_source,
        missing=len(need_lookup),
        unique=len(initial_need_lookup),
        attached=len(df) - len(need_lookup),
        uncovered=len(need_lookup),
        hierarchy_attached=int(hierarchy_attached),
    )

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
        df, attach_stats = attach_parent_molecule_ids(
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

    combined_stats = _merge_parent_stats(prep.parent_stats, attach_stats)

    logger.info(
        "parent_lookup_done",
        source=combined_stats.source,
        unique=combined_stats.unique,
        attached=combined_stats.attached,
        missing=combined_stats.missing,
        uncovered=combined_stats.uncovered,
        hierarchy_attached=combined_stats.hierarchy_attached,
        fallback_attached=combined_stats.fallback_attached,
        no_parent=combined_stats.no_parent,
    )

    return 0, ParentEnrichmentResult(df=df, parent_stats=combined_stats)


__all__ = [
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
]
