"""Command line interface for retrieving ChEMBL and PubChem compound data."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timedelta, timezone

# ruff: noqa: E402
from pathlib import Path

if __package__ is None:  # running as a script
    sys.path.append(str(Path(__file__).resolve().parents[1]))

import argparse
from collections import ChainMap
from collections.abc import Mapping, MutableMapping, Sequence
from dataclasses import dataclass
from itertools import islice
from typing import Any, Callable

import pandas as pd
import requests
from pandera.errors import SchemaErrors

from library import chembl_library as cl
from library import io
from library import molecule_catalog
from library import testitem_enrichment
from library import pubchem_library as pl
from library.molecule_catalog import (
    load_parent_catalog,
    query_parent_catalog,
    update_parent_catalog_cache,
    write_parent_catalog_cache,
)
from library.chembl_client import ChemblClient
from library.cli import (
    LoggerConfig,
    apply_config_overrides,
    configure_logger,
)
from library.cli import (
    build_parser as base_parser,
)
from library.config import (
    ApiCfg,
    Config,
    MoleculeCatalogCfg,
    PubChemCfg,
    _serialize_paths,
    ensure_dirs,
    print_config,
)
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.pipeline_metadata import add_pipeline_metadata
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from library.validation import validate_testitems
from schemas import TestitemsSchema, normalize_testitems

_TYPO_PARENT_COLUMN = "parant_molecule_id"


def ensure_no_parant_column(df: pd.DataFrame) -> None:
    """Raise a :class:`ValueError` if the legacy typo column is present."""

    if _TYPO_PARENT_COLUMN in df.columns:
        raise ValueError(
            "unexpected column 'parant_molecule_id'; use 'parent_molecule_id' instead"
        )


PARENT_LOOKUP_SOURCE_CACHE = "cache"
PARENT_LOOKUP_SOURCE_PARTIAL = "partial"
PARENT_LOOKUP_SOURCE_SYNC = "sync"
PARENT_LOOKUP_SOURCE_SKIPPED = "skipped"


PUBCHEM_COLUMNS = [
    "pubchem_cid",
    "pubchem_iupac_name",
    "pubchem_molecular_formula",
    "pubchem_isomeric_smiles",
    "pubchem_canonical_smiles",
    "pubchem_inchi",
    "pubchem_inchikey",
]


_PUBCHEM_CID_CACHE_ENCODING = "utf-8"
_CID_CACHE_MISSING = object()


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


_PUBCHEM_CACHE_SCHEMA_VERSION = 1


def _load_pubchem_cid_cache(
    path: Path | None, *, ttl_hours: float | None = None
) -> dict[str, str | None]:
    """Load CID cache mapping from disk."""

    if path is None:
        return {}
    try:
        with path.open("r", encoding=_PUBCHEM_CID_CACHE_ENCODING) as handle:
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
                updated_at = updated_at.replace(tzinfo=timezone.utc)
            if updated_at is not None:
                expires_at = updated_at + timedelta(hours=ttl_hours)
                if expires_at <= datetime.now(timezone.utc):
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


def _write_pubchem_cid_cache(path: Path | None, cache: Mapping[str, str | None]) -> None:
    """Persist CID cache mapping to disk."""

    if path is None:
        return
    ensure_dirs(path.parent)
    serialisable: dict[str, str] = {}
    for key, value in cache.items():
        if not key:
            continue
        if value is None:
            continue
        serialisable[key] = value
    try:
        with path.open("w", encoding=_PUBCHEM_CID_CACHE_ENCODING) as handle:
            payload = {
                "metadata": {
                    "version": _PUBCHEM_CACHE_SCHEMA_VERSION,
                    "updated_at": datetime.now(timezone.utc).isoformat(),
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


def _resolve_cid_from_row(
    row: pd.Series,
    cache: MutableMapping[str, str | None],
    cfg: PubChemCfg,
    *,
    chembl_id: str | None,
) -> str | None:
    """Return a CID for ``row`` using its structure fields."""


    if chembl_id and chembl_id in cache and cache[chembl_id]:
        return cache[chembl_id]
 


    chembl_id = _normalise_identifier(row.get("molecule_chembl_id"), uppercase=True)
 

    if chembl_id:
        cached = cache.get(chembl_id, _CID_CACHE_MISSING)
        if cached is not _CID_CACHE_MISSING:
            cid = _select_primary_cid(
                cached,
                chembl_id=chembl_id,
                identifier="cache",
                value=cached,
            )
            if chembl_id not in cache or cid != cached:
                cache[chembl_id] = cid
            return cid


    def _store(value: str | None) -> str | None:
        if chembl_id:
            cache[chembl_id] = value
        return value

    def _attempt(
        identifier: str,
        value: str | None,
        resolver: Callable[[str, PubChemCfg], str | None],
    ) -> str | None:
        if not value:
            return None
        resolved = resolver(value, cfg)
        return _select_primary_cid(
            resolved,
            chembl_id=chembl_id,
            identifier=identifier,
            value=value,
        )

    inchikey = _normalise_identifier(row.get("standard_inchi_key"), uppercase=True)
    cid = _attempt("standard_inchi_key", inchikey, pl.get_cid_from_inchikey)
    if cid:
        return _store(cid)

    inchi = _normalise_identifier(row.get("standard_inchi"))
    cid = _attempt("standard_inchi", inchi, pl.get_cid_from_inchi)
    if cid:
        return _store(cid)

    pref_name = _normalise_identifier(row.get("pref_name"))
    cid = _attempt("pref_name", pref_name, pl.get_cid)
    if cid:
        return _store(cid)
    cid = _attempt("pref_name_partial", pref_name, pl.get_all_cid)
    if cid:
        return _store(cid)

    smiles = _normalise_identifier(row.get("canonical_smiles"))
    cid = _attempt("canonical_smiles", smiles, pl.get_cid_from_smiles)
    if cid:
        return _store(cid)

    return None


def resolve_pubchem_cid(
    row: pd.Series,
    cache: MutableMapping[str, str | None],
    cfg: PubChemCfg,
    *,
    parent_loader: Callable[[str], pd.Series | None] | None = None,
) -> str | None:
    """Resolve PubChem CID for a ChEMBL record."""

    chembl_id = _normalise_identifier(row.get("molecule_chembl_id"), uppercase=True)
    cid = _resolve_cid_from_row(row, cache, cfg, chembl_id=chembl_id)
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

    parent_cid = _resolve_cid_from_row(parent_row, cache, cfg, chembl_id=parent_id)
    if parent_cid:
        cache[parent_id] = parent_cid
        if chembl_id:
            cache[chembl_id] = parent_cid
        return parent_cid


    return None


def resolve_pubchem_cid(
    row: pd.Series,
    cache: MutableMapping[str, str | None],
    cfg: PubChemCfg,
    *,
    parent_loader: Callable[[str], pd.Series | None] | None = None,
) -> str | None:
    """Resolve PubChem CID for a ChEMBL record."""

    chembl_id = _normalise_identifier(row.get("molecule_chembl_id"), uppercase=True)
    cid = _resolve_cid_from_row(row, cache, cfg, chembl_id=chembl_id)
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

    parent_cid = _resolve_cid_from_row(parent_row, cache, cfg, chembl_id=parent_id)
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

    normalised = (
        series.fillna("")
        .astype("string")
        .str.strip()
        .str.upper()
    )
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
) -> tuple[pd.DataFrame, ParentLookupStats]:
    """Attach parent molecule identifiers using the ChEMBL catalogue.

    The function operates on the original ``df`` content and can enrich frames
    that already contain parent identifiers as well as frames without the
    ``catalog_cfg.parent_field`` column.
    """

    if "parant_molecule_id" in df.columns:
        raise ValueError(
            "Input frame contains unexpected column 'parant_molecule_id'."
        )

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
        )
        return result, stats

    source_resolved = source
    normalised_child = _normalise_chembl_ids(result[child_column])
    unique_children = normalised_child[normalised_child != ""].unique()
    catalog_data: MutableMapping[str, str]
    used_partial_cache = False
    needs_full_sync = False
    partial_fetch_used = False
    full_sync_used = False

    if catalog is not None:
        base_view = {
            key: catalog[key]
            for key in unique_children
            if key in catalog
        }
        catalog_data = ChainMap({}, base_view)
    else:
        catalog_data = {}
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
            else:
                needs_full_sync = True

    parent_map = {
        key: catalog_data[key]
        for key in unique_children
        if key in catalog_data
    }
    missing_ids = [key for key in unique_children if key not in parent_map]

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
            if used_partial_cache:
                update_parent_catalog_cache(fetched, catalog_cfg)
            else:
                write_parent_catalog_cache(catalog_data, catalog_cfg)
                used_partial_cache = True

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

    parent_series = normalised_child.map(parent_map).astype("string")

    if parent_column in result.columns:
        existing_parent = result[parent_column].astype("string")
    else:
        existing_parent = pd.Series(pd.NA, index=result.index, dtype="string")

    combined_parent = existing_parent.copy()
    update_mask = combined_parent.isna() | combined_parent.eq("")
    combined_parent.loc[update_mask] = parent_series.loc[update_mask]
    result[parent_column] = combined_parent.astype("string")

    missing = int(combined_parent.isna().sum())
    attached = len(result) - missing

    final_source = source_resolved
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
    )

    logger.info(
        "parent_lookup_progress",
        source=stats.source,
        unique=stats.unique,
        attached=stats.attached,
        missing=stats.missing,
    )

    return result, stats


def add_pubchem_data(
    df: pd.DataFrame,
    cfg: PubChemCfg,
    *,
    client: ChemblClient | None = None,
    api_cfg: ApiCfg | None = None,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Augment ChEMBL records with PubChem information.

    For each canonical SMILES string in ``df``, the function looks up the
    corresponding PubChem CID and basic chemical properties. The PubChem
    fields are appended to the input frame. If a SMILES string cannot be
    resolved, empty values are inserted.

    Parameters
    ----------
    df:
        Data frame returned by :func:`library.chembl_library.get_testitem`.
    cfg:
        PubChem configuration options.
    client:
        Optional ChEMBL client reused to fetch parent structures.
    api_cfg:
        Optional API configuration required when ``client`` is supplied.
    timeout:
        Optional timeout for parent lookups via :func:`cl.get_testitem`.

    Returns
    -------
    pandas.DataFrame
        ``df`` with additional PubChem columns.

    """
    if df.empty:
        return df

    result = df.reset_index(drop=True).copy()
    prefer_local = getattr(cfg, "prefer_local_smiles", False)
    existing_cols = [col for col in PUBCHEM_COLUMNS if col in result.columns]
    if prefer_local and existing_cols:
        normalised = pd.DataFrame(
            {col: result[col].astype("string").replace("", pd.NA) for col in existing_cols}
        )
        complete_mask = normalised.notna().all(axis=1)
    else:
        complete_mask = pd.Series(False, index=result.index)

    cache_path = getattr(cfg, "cid_cache_path", None)
    cache_ttl_hours = getattr(cfg, "cache_ttl_hours", None)
    cid_cache = _load_pubchem_cid_cache(cache_path, ttl_hours=cache_ttl_hours)
    cache_dirty = False

    local_records: dict[str, pd.Series] = {}
    if "molecule_chembl_id" in result.columns:
        for index, value in result["molecule_chembl_id"].items():
            chembl_norm = _normalise_identifier(value, uppercase=True)
            if chembl_norm and chembl_norm not in local_records:
                local_records[chembl_norm] = result.loc[index]

    parent_record_cache: dict[str, pd.Series | None] = {}

    def load_parent_record(parent_id: str) -> pd.Series | None:
        parent_norm = _normalise_identifier(parent_id, uppercase=True)
        if not parent_norm:
            return None
        if parent_norm in parent_record_cache:
            return parent_record_cache[parent_norm]
        if parent_norm in local_records:
            parent_record_cache[parent_norm] = local_records[parent_norm]
            return parent_record_cache[parent_norm]
        if client is None or api_cfg is None:
            parent_record_cache[parent_norm] = None
            return None
        try:
            parent_df = cl.get_testitem(
                [parent_norm],
                cfg=api_cfg,
                client=client,
                chunk_size=1,
                timeout=timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.warning(
                "pubchem_parent_load_failed",
                parent=parent_norm,
                error=str(exc),
            )
            parent_record_cache[parent_norm] = None
            return None
        if parent_df.empty:
            parent_record_cache[parent_norm] = None
            return None
        parent_row = parent_df.iloc[0]
        local_records[parent_norm] = parent_row
        parent_record_cache[parent_norm] = parent_row
        return parent_row

    if "molecule_chembl_id" in result.columns and "pubchem_cid" in result.columns:
        for chembl_raw, cid_raw in zip(
            result["molecule_chembl_id"], result["pubchem_cid"]
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

    total = sum(
        1
        for idx, row in result.iterrows()
        if not (
            (chembl_id := _normalise_identifier(row.get("molecule_chembl_id"), uppercase=True))
            and chembl_id in cid_cache
            and cid_cache[chembl_id]
        )
    )
    if total:
        logger.info("pubchem_start", total=total)
    else:
        logger.info("pubchem_no_smiles")

    cid_by_index: dict[int, str | None] = {}
    progress = 0
    for idx, row in result.iterrows():
        if prefer_local and bool(complete_mask.loc[idx]):
            cid_by_index[idx] = _normalise_identifier(row.get("pubchem_cid"))
            continue

        chembl_id = _normalise_identifier(row.get("molecule_chembl_id"), uppercase=True)
        lookup_required = not (
            chembl_id and chembl_id in cid_cache and cid_cache[chembl_id]
        )
        before_present = chembl_id in cid_cache if chembl_id else False
        before_value = cid_cache[chembl_id] if before_present else _CID_CACHE_MISSING
        if lookup_required and total:
            progress += 1
            logger.info("pubchem_progress", current=progress, total=total)
        cid = resolve_pubchem_cid(
            row,
            cid_cache,
            cfg,
            parent_loader=load_parent_record,
        )
        cid_by_index[idx] = cid
        if chembl_id:
            after_present = chembl_id in cid_cache
            after_value = cid_cache[chembl_id] if after_present else _CID_CACHE_MISSING
            if before_present != after_present or after_value != before_value:
                cache_dirty = True

    if cache_dirty:
        _write_pubchem_cid_cache(cache_path, cid_cache)

    cid_list = [cid_by_index.get(idx) for idx in result.index]

    def _value_or_na(value: str | None) -> object:
        return value if value not in (None, "") else pd.NA

    lookup_cids: set[str] = set()
    for idx, cid in zip(result.index, cid_list):
        if not cid:
            continue
        if prefer_local and bool(complete_mask.loc[idx]):
            continue
        lookup_cids.add(cid)

    properties: dict[str, pl.Properties] = {}
    for cid in sorted(lookup_cids):
        properties[cid] = pl.get_properties(cid, cfg)

    pubchem_rows: list[dict[str, object]] = []
    for cid in cid_list:
        row_data: dict[str, object] = {col: pd.NA for col in PUBCHEM_COLUMNS}
        if cid:
            row_data["pubchem_cid"] = cid
            props = properties.get(cid)
            if props:
                row_data["pubchem_iupac_name"] = _value_or_na(props.IUPACName)
                row_data["pubchem_molecular_formula"] = _value_or_na(
                    props.MolecularFormula
                )
                row_data["pubchem_isomeric_smiles"] = _value_or_na(props.iSMILES)
                row_data["pubchem_canonical_smiles"] = _value_or_na(props.cSMILES)
                row_data["pubchem_inchi"] = _value_or_na(props.InChI)
                row_data["pubchem_inchikey"] = _value_or_na(props.InChIKey)
        pubchem_rows.append(row_data)

    pubchem_df = pd.DataFrame(pubchem_rows, index=result.index).convert_dtypes()

    prefer_values = getattr(cfg, "prefer_local_values", False)

    for col in PUBCHEM_COLUMNS:
        if col not in pubchem_df.columns:
            continue
        new_series = pubchem_df[col].astype("string")
        if col in result.columns:
            existing = result[col].astype("string")
            result[col] = existing
            if prefer_local:
                missing_mask = existing.isna() | existing.eq("")
                result.loc[missing_mask, col] = new_series[missing_mask]
            elif prefer_values:
                incoming = new_series.where(~(new_series.isna() | new_series.eq("")), pd.NA)
                combined = incoming.combine_first(existing)
                result[col] = combined.astype("string")
            else:
                result[col] = new_series
        else:
            result[col] = new_series

    return result


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute compound retrieval from the ChEMBL API and augment with PubChem data.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    limit = cfg.testitem.limit
    if limit is not None and limit < 0:
        logger.error(
            "invalid_limit",
            section="testitem.limit",
            limit=limit,
        )
        return 1

    # Initialise HTTP sessions for downstream HTTP calls
    pl.init_session(cfg.api, cfg.retry)
    # Initialise HTTP session for subsequent ChEMBL requests
    parent_stats = ParentLookupStats(
        source=PARENT_LOOKUP_SOURCE_SKIPPED,
        missing=0,
        unique=0,
        attached=0,
    )
    with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
        try:
            ids_iter = io.read_ids(
                args.input_csv, column=cfg.testitem.column, cfg=cfg.io
            )
            if limit is not None:
                ids = list(islice(ids_iter, limit))
                logger.info("process_limit", limit=len(ids))
            else:
                ids = list(ids_iter)
        except (FileNotFoundError, ValueError) as exc:
            logger.error(
                "read_fail",
                error=str(exc),
                path=str(args.input_csv),
            )
            return 1

        logger.info("identifiers_retrieved", count=len(ids))
        logger.info("chembl_fetch_start", batch_size=cfg.testitem.batch_size)

        try:
            df = cl.get_testitem(
                ids,
                cfg=cfg.api,
                client=client,
                chunk_size=cfg.testitem.batch_size,
                timeout=cfg.testitem.timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error(
                "testitem_fetch_failed",
                error=str(exc),
                batch_size=cfg.testitem.batch_size,
                timeout=cfg.testitem.timeout,
            )
            return 1
        logger.info("chembl_fetch_done", rows=len(df))
        parent_column = cfg.molecule_catalog.parent_field
        child_column = cfg.molecule_catalog.child_field

        if child_column in df.columns:
            normalised_ids = _normalise_chembl_ids(df[child_column])
        else:
            normalised_ids = pd.Series("", index=df.index, dtype="string")

        if parent_column in df.columns:
            existing_parent = _normalise_chembl_ids(df[parent_column])
        else:
            existing_parent = pd.Series("", index=df.index, dtype="string")

        need_lookup_mask = (normalised_ids != "") & (existing_parent == "")
        need_lookup = set(normalised_ids[need_lookup_mask])

        cache_before = _cache_state(cfg.molecule_catalog.cache_path)
        cache_after = cache_before
        parent_catalog: dict[str, str] = {}
        parent_catalog_source = PARENT_LOOKUP_SOURCE_SKIPPED

        if need_lookup and cache_before[0]:
            try:
                parent_catalog = query_parent_catalog(
                    need_lookup,
                    catalog_cfg=cfg.molecule_catalog,
                )
            except (requests.RequestException, ValueError) as exc:
                logger.error(
                    "parent_catalog_invalid",
                    error=str(exc),
                    path=str(cfg.molecule_catalog.cache_path),
                )
                return 1
            cache_after = _cache_state(cfg.molecule_catalog.cache_path)
            parent_catalog_source = (
                PARENT_LOOKUP_SOURCE_CACHE
                if cache_after == cache_before
                else PARENT_LOOKUP_SOURCE_SYNC
            )
            if parent_catalog:
                need_lookup -= set(parent_catalog)

        fetched_remote = False
        if need_lookup and not parent_catalog:
            try:
                fallback_catalog = load_parent_catalog(
                    client=client,
                    api_cfg=cfg.api,
                    catalog_cfg=cfg.molecule_catalog,
                    timeout=cfg.testitem.timeout,
                )
            except (requests.RequestException, ValueError) as exc:
                logger.error("parent_catalog_invalid", error=str(exc))
                return 1
            if fallback_catalog:
                parent_catalog.update(fallback_catalog)
                parent_catalog_source = PARENT_LOOKUP_SOURCE_CACHE
                need_lookup -= set(fallback_catalog)
        if need_lookup:
            try:
                fetched = molecule_catalog.fetch_parent_catalog_for(
                    need_lookup,
                    client=client,
                    api_cfg=cfg.api,
                    timeout=cfg.testitem.timeout,
                    catalog_cfg=cfg.molecule_catalog,
                )
            except (requests.RequestException, ValueError) as exc:
                logger.error("parent_lookup_partial_fetch_failed", error=str(exc))
                return 1
            if fetched:
                parent_catalog.update(fetched)
                update_parent_catalog_cache(fetched, cfg.molecule_catalog)
                parent_catalog_source = PARENT_LOOKUP_SOURCE_PARTIAL
        try:
            ensure_no_parant_column(df)
        except ValueError as exc:
            logger.error(
                "invalid_column",
                column=_TYPO_PARENT_COLUMN,
                error=str(exc),
            )
            return 1
 
        logger.info("parent_lookup_start")
        try:
            catalog_arg = parent_catalog if parent_catalog else (
                None if need_lookup else parent_catalog
            )
            df, parent_stats = attach_parent_molecule_ids(
                df,
                client=client,
                api_cfg=cfg.api,
                catalog_cfg=cfg.molecule_catalog,
                timeout=cfg.testitem.timeout,
                catalog=catalog_arg,
                source=parent_catalog_source,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error("parent_lookup_failed", error=str(exc))
            return 1
        logger.info(
            "parent_lookup_done",
            source=parent_stats.source,
            unique=parent_stats.unique,
            attached=parent_stats.attached,
            missing=parent_stats.missing,
        )

        logger.info("pubchem_augment_start")
        df = add_pubchem_data(
            df,
            cfg.pubchem,
            client=client,
            api_cfg=cfg.api,
            timeout=cfg.testitem.timeout,
        )
        logger.info("pubchem_augment_done")

        enrichment_cfg = cfg.testitem_molecule_enrichment
        if enrichment_cfg.enable:
            logger.info("testitem_enrichment_start")
            try:
                df = testitem_enrichment.enrich(
                    df,
                    cfg=enrichment_cfg,
                    io_cfg=cfg.io,
                )
            except ValueError as exc:
                logger.error("testitem_enrichment_failed", error=str(exc))
                return 1
            logger.info("testitem_enrichment_done")

        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        df = normalize_testitems(df)
        df = add_pipeline_metadata(df)
        # Determine column order: schema columns first, followed by
        # additional fields sorted alphabetically.
        schema_cols = list(TestitemsSchema.columns)
        head = [c for c in schema_cols if c in df.columns]
        tail = sorted(c for c in df.columns if c not in schema_cols)
        col_order = head + tail
        rows_total = len(df)
        exit_code = 0
    required_cols = {
        name for name, col in TestitemsSchema.columns.items() if col.required
    }
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
        csv_path = io.write_csv(
            df,
            output,
            cfg=cfg,
            key_cols=key_cols or None,
            col_order=col_order,
        )
        logger.info("write_done", rows=rows_kept, path=str(csv_path))
    except OSError as exc:
        logger.error(
            "write_fail",
            error=str(exc),
            path=str(output),
        )
        return 1

    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_kept,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(csv_path),
        "parent_lookup_source": parent_stats.source,
        "parent_lookup_missing": parent_stats.missing,
    }
    write_meta_yaml(
        csv_path=csv_path,
        command=" ".join(sys.argv),
        config_subset=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(args.input_csv)},
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


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser, log_cfg = base_parser(
        "ChEMBL and PubChem compound data utilities",
        column="molecule_chembl_id",
        chunk_size=5,
        size_option="--batch-size",
        size_dest="batch_size",
    )
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
        help="Maximum number of identifiers to process",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be a positive integer")
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "timeout": "testitem.timeout",
                "column": "testitem.column",
                "batch_size": "testitem.batch_size",
                "limit": "testitem.limit",
            },
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error(
            "directory_setup_failed",
            error=str(exc),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    exit_code: int = args.func(cfg, args)
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
