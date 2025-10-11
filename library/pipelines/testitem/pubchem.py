"""PubChem enrichment helpers for the test item pipeline."""

from __future__ import annotations

import json
import threading
from collections.abc import (
    Callable,
    Collection,
    Hashable,
    Mapping,
    MutableMapping,
    Sequence,
)
from concurrent.futures import CancelledError, ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from functools import lru_cache
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    TypeAlias,
    cast,
)

import pandas as pd
import requests

from library.common.log import logger
from library.config import ApiCfg, PubChemCfg, RetryCfg
from library.integration.chembl_client import ChemblClient
from library.utils.atomic import open_atomic

if TYPE_CHECKING:
    from library.integration.pubchem_library import Properties, PubChemResolution

UTC = timezone.utc  # noqa: UP017
ENCODING_UTF8 = "utf-8"
PUBCHEM_CID_CACHE_ENCODING = ENCODING_UTF8
_TESTITEM_PUBCHEM_COLUMNS = (
    "pubchem_cid",
    "pubchem_iupac_name",
    "pubchem_molecular_formula",
    "pubchem_isomeric_smiles",
    "pubchem_canonical_smiles",
    "pubchem_inchi",
    "pubchem_inchikey",
)
PUBCHEM_COLUMNS = list(_TESTITEM_PUBCHEM_COLUMNS)

_CID_CACHE_MISSING = object()
_PUBCHEM_CACHE_SCHEMA_VERSION = 1

_PUBCHEM_SESSION_SIGNATURE: str | None = None


ResolutionCache: TypeAlias = MutableMapping[Hashable, "PubChemResolution"]
_PUBCHEM_SESSION_LOCK = threading.Lock()


@lru_cache(maxsize=1)
def _load_chembl_library():
    """Return the ChemBL integration module without triggering circular imports."""

    from library.integration import chembl_library as chembl_lib

    return chembl_lib


@lru_cache(maxsize=1)
def _load_pubchem_library():
    """Return the PubChem integration module without triggering circular imports."""

    from library.integration import pubchem_library as pubchem_lib

    return pubchem_lib


class _PubChemLibraryProxy:
    """Proxy exposing ``library.integration.pubchem_library`` lazily."""

    def __getattr__(self, name: str) -> Any:
        return getattr(_load_pubchem_library(), name)

    def __setattr__(self, name: str, value: Any) -> None:
        setattr(_load_pubchem_library(), name, value)

    def __dir__(self) -> list[str]:  # pragma: no cover - debug helper
        return sorted(set(dir(type(self)) + dir(_load_pubchem_library())))


pl = _PubChemLibraryProxy()


def _pubchem_session_signature(api_cfg: ApiCfg, retry_cfg: RetryCfg) -> str:
    """Return a stable signature for the PubChem session configuration."""

    payload = {
        "api": api_cfg.model_dump(mode="json"),
        "retry": retry_cfg.model_dump(mode="json"),
    }
    return json.dumps(payload, sort_keys=True)


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


def _load_pubchem_cid_cache(
    path: Path | str | None, *, ttl_hours: float | None = None
) -> dict[str, str | None]:
    """Load CID cache mapping from disk."""

    if path is None:
        return {}
    cache_path = Path(path)
    try:
        with cache_path.open("r", encoding=PUBCHEM_CID_CACHE_ENCODING) as handle:
            data = json.load(handle)
    except FileNotFoundError:
        return {}
    except (OSError, json.JSONDecodeError) as exc:  # pragma: no cover - I/O errors
        logger.warning(
            "pubchem_cache_load_failed", path=str(cache_path), error=str(exc)
        )
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
        logger.warning("pubchem_cache_invalid", path=str(cache_path))
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
                    logger.info("pubchem_cache_expired", path=str(cache_path))
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
    path: Path | str | None, cache: Mapping[str, str | None]
) -> None:
    """Persist CID cache mapping to disk."""

    if path is None:
        return
    cache_path = Path(path)
    serialisable: dict[str, str | None] = {}
    for key, value in cache.items():
        if not key:
            continue
        if value is None:
            serialisable[key] = None
        else:
            serialisable[key] = value
    try:
        with open_atomic(cache_path, encoding=PUBCHEM_CID_CACHE_ENCODING) as handle:
            payload = {
                "metadata": {
                    "version": _PUBCHEM_CACHE_SCHEMA_VERSION,
                    "updated_at": datetime.now(UTC).isoformat(),
                },
                "values": serialisable,
            }
            json.dump(payload, handle, indent=2, sort_keys=True)
    except OSError as exc:  # pragma: no cover - I/O errors
        logger.warning(
            "pubchem_cache_write_failed", path=str(cache_path), error=str(exc)
        )


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
    resolution_cache: ResolutionCache | None = None,
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

    pubchem_lib = _load_pubchem_library()

    resolution = pubchem_lib.resolve_pubchem_record(
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

    temporary_failure = bool(getattr(resolution, "temporary_failure", False))

    if temporary_failure:
        return None

    if not cfg.use_parent_for_salts:
        if chembl_id and chembl_id not in cache and not temporary_failure:
            cache[chembl_id] = None
        return None

    parent_raw = None
    if isinstance(row, pd.Series) and "parent_molecule_chembl_id" in row.index:
        parent_raw = row.get("parent_molecule_chembl_id")
    parent_id = _normalise_identifier(parent_raw, uppercase=True)
    if not parent_id:
        if chembl_id and chembl_id not in cache and not temporary_failure:
            cache[chembl_id] = None
        return None

    if parent_loader is None:
        if chembl_id and chembl_id not in cache and not temporary_failure:
            cache[chembl_id] = None
        return None

    if parent_id in visited:
        logger.info(
            "pubchem_parent_structure_missing",
            child=chembl_id,
            parent=parent_id,
            reason="parent_cycle_detected",
        )
        if parent_id not in cache and not temporary_failure:
            cache[parent_id] = None
        if chembl_id and chembl_id not in cache and not temporary_failure:
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
        if chembl_id and chembl_id not in cache and not temporary_failure:
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
    if not temporary_failure:
        cache.setdefault(parent_id, None)
        cache.setdefault(chembl_id or parent_id, None)
    return None


def _prepare_pubchem_caches(
    frame: pd.DataFrame,
    cfg: PubChemCfg,
    *,
    cache_path: Path | str | None,
    cache_ttl_hours: float | None,
    cid_cache: MutableMapping[str, str | None] | None,
    resolution_cache: ResolutionCache | None,
    parent_record_cache: MutableMapping[str, pd.Series | None] | None,
) -> tuple[
    MutableMapping[str, str | None],
    ResolutionCache,
    MutableMapping[str, pd.Series | None],
    list[str],
    Callable[[str], pd.Series | None],
]:
    """Prepare caches and parent lookups required for PubChem enrichment."""

    if cid_cache is None:
        cid_cache = _load_pubchem_cid_cache(cache_path, ttl_hours=cache_ttl_hours)
    if resolution_cache is None:
        resolution_cache = {}
    if parent_record_cache is None:
        parent_record_cache = {}

    local_records: pd.DataFrame | None = None
    if "molecule_chembl_id" in frame.columns:
        prepared_local_records = frame.assign(
            __local_molecule=lambda data: data["molecule_chembl_id"]
            .astype("string")
            .str.strip()
            .str.upper()
        ).dropna(subset=["__local_molecule"])
        if not prepared_local_records.empty:
            local_records = (
                prepared_local_records.drop_duplicates("__local_molecule")
                .set_index("__local_molecule")
                .rename_axis("molecule_chembl_id")
                .reindex(columns=frame.columns)
            )
            for chembl_norm in local_records.index:
                try:
                    record = local_records.loc[chembl_norm]
                except KeyError:  # pragma: no cover - defensive
                    continue
                parent_record_cache[chembl_norm] = record.copy()

    pending_parent_ids: list[str] = []
    seen_parent_ids: set[str] = set()
    if "parent_molecule_chembl_id" in frame.columns:
        for parent_value in frame["parent_molecule_chembl_id"]:
            parent_norm = _normalise_identifier(parent_value, uppercase=True)
            if not parent_norm:
                continue
            if parent_norm in parent_record_cache:
                continue
            if local_records is not None and parent_norm in local_records.index:
                continue
            if parent_norm in seen_parent_ids:
                continue
            seen_parent_ids.add(parent_norm)
            pending_parent_ids.append(parent_norm)

    def load_parent_record(parent_norm: str) -> pd.Series | None:
        if parent_norm in parent_record_cache:
            return parent_record_cache[parent_norm]
        if local_records is not None and parent_norm in local_records.index:
            try:
                return local_records.loc[parent_norm]
            except KeyError:  # pragma: no cover - defensive
                return None
        return parent_record_cache.get(parent_norm)

    return (
        cid_cache,
        resolution_cache,
        parent_record_cache,
        pending_parent_ids,
        load_parent_record,
    )


def _prefetch_parents(
    parent_ids: Sequence[str],
    *,
    client: ChemblClient | None,
    api_cfg: ApiCfg | None,
    cfg: PubChemCfg,
    timeout: float | None,
    testitem_fields: Sequence[str] | None,
    request_limit: int,
    parent_record_cache: MutableMapping[str, pd.Series | None],
) -> None:
    """Prefetch parent records for the provided identifiers."""

    if not parent_ids or client is None or api_cfg is None:
        return

    logger.info(
        "pubchem_parent_prefetch",
        count=len(parent_ids),
        batch_size=getattr(cfg, "batch_size", 0),
    )
    chembl_lib = _load_chembl_library()

    try:
        fetched = chembl_lib.get_testitem(
            parent_ids,
            cfg=api_cfg,
            client=client,
            chunk_size=getattr(cfg, "batch_size", 0),
            timeout=timeout,
            fields=testitem_fields,
            page_limit=request_limit,
        )
    except (requests.RequestException, ValueError) as exc:
        logger.warning("pubchem_parent_prefetch_failed", error=str(exc))
        return

    if fetched.empty:
        return

    fetched = fetched.astype("string")
    fetched["molecule_chembl_id"] = fetched["molecule_chembl_id"].str.upper()
    fetched = fetched.drop_duplicates("molecule_chembl_id")
    fetched = fetched.set_index("molecule_chembl_id")
    for parent_norm, row in fetched.iterrows():
        parent_record_cache[parent_norm] = row


def _resolve_pubchem_cids(
    frame: pd.DataFrame,
    cfg: PubChemCfg,
    *,
    cid_cache: MutableMapping[str, str | None],
    resolution_cache: ResolutionCache,
    load_parent_record: Callable[[str], pd.Series | None] | None,
    skip_mask: pd.Series,
    prefer_local_mask: pd.Series,
) -> tuple[pd.Series, set[str], bool]:
    """Resolve PubChem CIDs for rows requiring lookup."""

    if "pubchem_cid" in frame.columns:
        cid_series = frame["pubchem_cid"].astype("string").copy()
    else:
        cid_series = pd.Series(pd.NA, index=frame.index, dtype="string")

    if "molecule_chembl_id" in frame.columns:
        chembl_norm = frame["molecule_chembl_id"].map(
            lambda value: _normalise_identifier(value, uppercase=True)
        )
    else:
        chembl_norm = pd.Series([None] * len(frame), index=frame.index, dtype="object")

    def _is_cached(chembl_id: str | None) -> bool:
        return bool(chembl_id) and chembl_id in cid_cache

    cached_mask = chembl_norm.map(_is_cached)
    needs_lookup_mask = (
        chembl_norm.notna()
        & ~skip_mask.astype(bool)
        & ~prefer_local_mask.astype(bool)
        & ~cached_mask
    )

    total = int(needs_lookup_mask.sum())
    if total:
        logger.info("pubchem_start", total=total)
    else:
        logger.info("pubchem_no_smiles")

    cid_series = cid_series.astype("string")
    lookup_cids: set[str] = set()
    cache_dirty = False
    for idx, chembl_id in chembl_norm[cached_mask].items():
        cached_value = cid_cache.get(chembl_id)
        if cached_value is None:
            continue
        cid_series.loc[idx] = cached_value
        lookup_cids.add(cached_value)

    for progress, row in enumerate(frame.loc[needs_lookup_mask].itertuples(), start=1):
        logger.info("pubchem_progress", current=progress, total=total)
        idx = row.Index
        chembl_id = chembl_norm.loc[idx]
        before_present = bool(chembl_id and chembl_id in cid_cache)
        before_value = cid_cache[chembl_id] if before_present else _CID_CACHE_MISSING
        cid = resolve_pubchem_cid(
            frame.loc[idx],
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

    return cid_series, lookup_cids, cache_dirty


def _merge_pubchem_properties(
    frame: pd.DataFrame,
    cid_series: pd.Series,
    lookup_cids: Collection[str],
    *,
    cfg: PubChemCfg,
    skip_mask: pd.Series,
    prefer_local_mask: pd.Series,
) -> pd.DataFrame:
    """Fetch PubChem properties and merge them with the provided frame."""

    def _value_or_na(value: str | None) -> object:
        return value if value not in (None, "") else pd.NA

    properties_records: dict[str, dict[str, object]] = {}

    lookup_order = sorted(lookup_cids)
    if lookup_order:
        configured_batch_size = max(int(getattr(cfg, "batch_size", 1)), 1)
        rps_limit = int(getattr(cfg, "rps", configured_batch_size))
        max_workers = max(1, min(configured_batch_size, rps_limit))
        batch_size = max_workers

        pubchem_lib = _load_pubchem_library()

        def _empty_properties() -> Properties:
            return pubchem_lib.Properties(None, None, None, None, None, None)

        def _fetch_properties(cid: str) -> Properties:
            return pubchem_lib.get_properties(cid, cfg)

        service_unavailable = False
        unavailable_error: str | None = None
        failed_cids: set[str] = set()

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for start in range(0, len(lookup_order), batch_size):
                if service_unavailable:
                    break
                batch = lookup_order[start : start + batch_size]
                future_map = {
                    executor.submit(_fetch_properties, cid): cid for cid in batch
                }
                for future, cid in future_map.items():
                    if service_unavailable:
                        cancelled = future.cancel()
                        if not cancelled:
                            try:
                                future.result()
                            except (  # pragma: no cover - defensive cleanup
                                CancelledError,
                                pubchem_lib.PubChemServiceUnavailable,
                                Exception,
                            ):
                                pass
                        props = _empty_properties()
                        failed_cids.add(cid)
                    else:
                        try:
                            props = future.result()
                        except pubchem_lib.PubChemServiceUnavailable as exc:
                            logger.warning(
                                "pubchem_properties_failed", cid=cid, error=str(exc)
                            )
                            service_unavailable = True
                            unavailable_error = str(exc)
                            failed_cids.add(cid)
                            props = _empty_properties()
                        except Exception as exc:  # pragma: no cover - defensive
                            logger.warning(
                                "pubchem_properties_failed", cid=cid, error=str(exc)
                            )
                            props = _empty_properties()
                    values = {
                        "pubchem_iupac_name": _value_or_na(props.IUPACName),
                        "pubchem_molecular_formula": _value_or_na(
                            props.MolecularFormula
                        ),
                        "pubchem_isomeric_smiles": _value_or_na(props.iSMILES),
                        "pubchem_canonical_smiles": _value_or_na(props.cSMILES),
                        "pubchem_inchi": _value_or_na(props.InChI),
                        "pubchem_inchikey": _value_or_na(props.InChIKey),
                    }
                    properties_records[cid] = values

                if service_unavailable:
                    remaining_cids = lookup_order[start + batch_size :]
                    if remaining_cids:
                        failed_cids.update(remaining_cids)
                    break

        if service_unavailable:
            pending = max(len(failed_cids), 0)
            log_payload = {"pending": pending}
            if unavailable_error:
                log_payload["error"] = unavailable_error
            logger.warning("pubchem_properties_unavailable", **log_payload)

    properties_df = pd.DataFrame.from_dict(properties_records, orient="index")
    pubchem_df = cid_series.to_frame("pubchem_cid").join(
        properties_df, on="pubchem_cid"
    )
    pubchem_df = pubchem_df.reindex(frame.index)

    preserve_mask = skip_mask | prefer_local_mask
    existing_columns = [column for column in PUBCHEM_COLUMNS if column in frame.columns]
    prefer_local_values = bool(getattr(cfg, "prefer_local_values", False))

    if existing_columns:
        original_existing = frame[existing_columns].astype("string")
        for column in existing_columns:
            if column not in pubchem_df.columns:
                pubchem_df[column] = pd.Series(
                    pd.NA, index=pubchem_df.index, dtype="string"
                )

        updated = pubchem_df[existing_columns].astype("string")

        if prefer_local_values:
            for column in existing_columns:
                new_col = updated[column]
                original_col = original_existing[column]
                new_missing = new_col.isna() | (new_col.fillna("").str.len() == 0)
                original_present = ~original_col.isna() & (
                    original_col.fillna("").str.len() > 0
                )
                replace_mask = new_missing & original_present
                if replace_mask.any():
                    updated[column] = new_col.mask(replace_mask, original_col)

        if preserve_mask.any():
            for column in existing_columns:
                updated[column] = updated[column].mask(
                    preserve_mask, original_existing[column]
                )

        pubchem_df[existing_columns] = updated

    return pubchem_df.convert_dtypes()


def add_pubchem_data(
    df: pd.DataFrame,
    cfg: PubChemCfg,
    *,
    client: ChemblClient | None = None,
    api_cfg: ApiCfg | None = None,
    timeout: float | None = None,
    cid_cache: MutableMapping[str, str | None] | None = None,
    resolution_cache: ResolutionCache | None = None,
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
    (
        cid_cache,
        resolution_cache,
        parent_record_cache,
        pending_parent_ids,
        load_parent_record,
    ) = _prepare_pubchem_caches(
        result,
        cfg,
        cache_path=cache_path,
        cache_ttl_hours=cache_ttl_hours,
        cid_cache=cid_cache,
        resolution_cache=resolution_cache,
        parent_record_cache=parent_record_cache,
    )
    cache_dirty = False

    _prefetch_parents(
        pending_parent_ids,
        client=client,
        api_cfg=api_cfg,
        cfg=cfg,
        timeout=timeout,
        testitem_fields=testitem_fields,
        request_limit=request_limit,
        parent_record_cache=parent_record_cache,
    )

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

    skip_mask = pd.Series(index=result.index, data=False, dtype="bool")
    if skip_indexes:
        skip_mask.loc[list(skip_indexes)] = True

    prefer_local_mask = complete_mask.astype(bool)

    cid_series, lookup_cids, cache_modified = _resolve_pubchem_cids(
        result,
        cfg,
        cid_cache=cid_cache,
        resolution_cache=resolution_cache,
        load_parent_record=load_parent_record,
        skip_mask=skip_mask,
        prefer_local_mask=prefer_local_mask,
    )
    cache_dirty = cache_dirty or cache_modified

    pubchem_df = _merge_pubchem_properties(
        result,
        cid_series,
        lookup_cids,
        cfg=cfg,
        skip_mask=skip_mask,
        prefer_local_mask=prefer_local_mask,
    )

    result = result.drop(
        columns=[col for col in PUBCHEM_COLUMNS if col in result.columns]
    )
    result = result.join(pubchem_df)

    if cache_dirty:
        _write_pubchem_cid_cache(cache_path, cid_cache)

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
    """Augment ``df`` with PubChem information if enabled."""

    pubchem_cid_cache: dict[str, str | None] | None = None
    pubchem_resolution_cache: ResolutionCache | None = None
    pubchem_parent_record_cache: dict[str, pd.Series | None] | None = None
    if getattr(pubchem_cfg, "enable", True):
        global _PUBCHEM_SESSION_SIGNATURE
        session_signature = _pubchem_session_signature(api_cfg, retry_cfg)
        with _PUBCHEM_SESSION_LOCK:
            if session_signature != _PUBCHEM_SESSION_SIGNATURE:
                _load_pubchem_library().init_session(api_cfg, retry_cfg)
                _PUBCHEM_SESSION_SIGNATURE = session_signature
        pubchem_cid_cache = _load_pubchem_cid_cache(
            getattr(pubchem_cfg, "cid_cache_path", None),
            ttl_hours=getattr(pubchem_cfg, "cache_ttl_hours", None),
        )
        pubchem_resolution_cache = cast(ResolutionCache, {})
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


__all__ = [
    "PUBCHEM_CID_CACHE_ENCODING",
    "PUBCHEM_COLUMNS",
    "_CID_CACHE_MISSING",
    "_PUBCHEM_CACHE_SCHEMA_VERSION",
    "_PUBCHEM_SESSION_LOCK",
    "_PUBCHEM_SESSION_SIGNATURE",
    "_prepare_pubchem_caches",
    "_prefetch_parents",
    "_pubchem_session_signature",
    "_resolve_pubchem_cids",
    "_merge_pubchem_properties",
    "_write_pubchem_cid_cache",
    "_load_pubchem_cid_cache",
    "add_pubchem_data",
    "augment_pubchem",
    "resolve_pubchem_cid",
    "pl",
]
