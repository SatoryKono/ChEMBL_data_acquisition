"""Utilities for retrieving and caching the ChEMBL molecule parent catalogue."""

from __future__ import annotations

import csv
import json
import sqlite3
import threading
from collections.abc import Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from urllib.parse import urlencode, urljoin

from library.clients import ChemblClient, _chunked
from library.utils.atomic import open_atomic

from ..common.log import logger
from ..common.rate_limiter import sleep
from ..config import ApiCfg, MoleculeCatalogCfg

_DEFAULT_CATALOG_CFG = MoleculeCatalogCfg()
_PARENT_LOOKUP_FALLBACK_THRESHOLD = 1
_PARENT_LOOKUP_RETRY_THRESHOLD = 10
_PARENT_LOOKUP_SINGLE_ATTEMPTS_MIN = 1
_SQLITE_VARIABLE_LIMIT = 900
_PARENT_LOOKUP_SINGLE_CONCURRENCY = 4

_PARENT_CATALOG_LOCK = threading.Lock()

__all__ = [
    "fetch_parent_catalog",
    "fetch_parent_catalog_for",
    "fetch_parent_for_id",
    "load_parent_catalog",
    "query_parent_catalog",
    "update_parent_catalog_cache",
    "write_parent_catalog_cache",
]


def _resolve_lookup_fields(catalog_cfg: MoleculeCatalogCfg) -> tuple[str, ...]:
    fields = tuple(catalog_cfg.fields or ())
    if fields:
        return fields
    return (catalog_cfg.child_field, catalog_cfg.parent_field)


def _normalise_chembl_id(value: str) -> str:
    """Return *value* normalised as an upper-case ChEMBL identifier."""

    normalised = value.strip().upper()
    if not normalised:
        raise ValueError("empty identifier")
    return normalised


def _filters_exclude_parentless(catalog_cfg: MoleculeCatalogCfg) -> bool:
    flag = catalog_cfg.filters.get(f"{catalog_cfg.parent_field}__isnull")
    if flag is None:
        return False
    if isinstance(flag, bool):
        return flag is False
    return str(flag).strip().lower() in {"false", "0"}


def _build_initial_url(api_cfg: ApiCfg, catalog_cfg: MoleculeCatalogCfg) -> str:
    base = api_cfg.chembl_base.rstrip("/")
    resource = catalog_cfg.endpoint.lstrip("/")
    if not resource.endswith(".json"):
        resource = f"{resource}.json"
    params: dict[str, str] = {
        "format": "json",
        "limit": str(catalog_cfg.page_size),
    }
    if catalog_cfg.fields:
        params["fields"] = ",".join(catalog_cfg.fields)
    for key, value in catalog_cfg.filters.items():
        params[key] = str(value)
    query = urlencode(params)
    return f"{base}/{resource}?{query}"


def fetch_parent_for_id(
    chembl_id: str,
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    timeout: float | None = None,
    catalog_cfg: MoleculeCatalogCfg | None = None,
) -> tuple[str, str] | None:
    """Return a single child-to-parent mapping for ``chembl_id``.

    Parameters
    ----------
    chembl_id:
        Identifier of the molecule to query.
    client:
        HTTP client used to perform the request.
    api_cfg:
        API configuration used for timeouts and retry settings.
    timeout:
        Optional override for the read timeout.

    Returns
    -------
    tuple[str, str] | None
        A tuple containing the normalised child and parent identifiers when
        available; otherwise ``None``.
    """

    try:
        child_id = _normalise_chembl_id(chembl_id)
    except ValueError:
        return None

    cfg = catalog_cfg or _DEFAULT_CATALOG_CFG
    base = api_cfg.chembl_base.rstrip("/")
    fields = _resolve_lookup_fields(cfg)
    params = {
        "format": "json",
        "fields": ",".join(fields),
    }
    url = f"{base}/molecule/{child_id}.json?{urlencode(params)}"
    effective_timeout = timeout if timeout is not None else api_cfg.timeout_read
    data = client.request_json(url, cfg=api_cfg, timeout=effective_timeout)

    item = (
        data.get("molecule")
        or data.get("molecule_parent")
        or data.get("molecule_parents")
        or data
    )
    if isinstance(item, list):
        item = item[0] if item else None

    if not isinstance(item, Mapping):
        logger.warning("unexpected_single_response", extra={"url": url})
        return None

    child_value = item.get(cfg.child_field) or child_id
    parent_value = item.get(cfg.parent_field)
    if not parent_value:
        return None

    try:
        child_norm = _normalise_chembl_id(str(child_value))
        parent_norm = _normalise_chembl_id(str(parent_value))
    except ValueError:
        return None

    return child_norm, parent_norm


def _fetch_parent_catalog_chunk(
    chunk: list[str],
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    timeout: float | None,
    catalog_cfg: MoleculeCatalogCfg,
) -> dict[str, str]:
    fields = _resolve_lookup_fields(catalog_cfg)
    params = {
        "format": "json",
        "limit": str(len(chunk)),
        f"{catalog_cfg.child_field}__in": ",".join(chunk),
        "fields": ",".join(fields),
    }
    for key, value in catalog_cfg.filters.items():
        params[key] = str(value)
    base = api_cfg.chembl_base.rstrip("/")
    url = f"{base}/molecule.json?{urlencode(params)}"
    effective_timeout = timeout if timeout is not None else api_cfg.timeout_read
    data = client.request_json(url, cfg=api_cfg, timeout=effective_timeout)
    items_raw = (
        data.get("molecules")
        or data.get("molecule")
        or data.get("molecule_parents")
        or data.get("molecule_parent")
        or []
    )
    if isinstance(items_raw, Mapping):
        items = [items_raw]
    elif isinstance(items_raw, list):
        items = items_raw
    else:
        logger.warning("unexpected_response_items", extra={"url": url})
        return {}

    allowed = set(chunk)
    result: dict[str, str] = {}
    for item in items:
        if not isinstance(item, Mapping):
            continue
        child = item.get(catalog_cfg.child_field)
        parent = item.get(catalog_cfg.parent_field)
        if not child or not parent:
            continue
        try:
            child_id = _normalise_chembl_id(str(child))
            parent_id = _normalise_chembl_id(str(parent))
        except ValueError:
            continue
        if child_id not in allowed:
            continue
        result[child_id] = parent_id

    return result


def _fetch_parent_catalog_via_helper(
    ids: Iterable[str],
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    timeout: float | None,
    existing: Mapping[str, str],
    catalog_cfg: MoleculeCatalogCfg | None = None,
    allow_rebatch: bool = True,
    allow_single_lookup: bool = True,
) -> dict[str, str]:
    cfg = catalog_cfg or _DEFAULT_CATALOG_CFG
    pending: list[str] = []
    seen: set[str] = set()
    for chembl_id in ids:
        if chembl_id in existing or chembl_id in seen:
            continue
        seen.add(chembl_id)
        pending.append(chembl_id)

    if not pending:
        return {}

    attempts = max(_PARENT_LOOKUP_SINGLE_ATTEMPTS_MIN, api_cfg.retries + 1)
    delay = api_cfg.backoff_factor
    single_limit = cfg.fallback_single_limit
    result: dict[str, str] = {}
    batch_attempts = 0
    single_requests = 0
    fallback_start = perf_counter()

    try:
        if not allow_single_lookup:
            if allow_rebatch:
                chunk_size = max(1, cfg.page_size)
            else:
                chunk_size = max(1, len(pending))
            for chunk in _chunked(pending, chunk_size):
                batch_attempts += 1
                try:
                    chunk_result = _fetch_parent_catalog_chunk(
                        chunk,
                        client=client,
                        api_cfg=api_cfg,
                        timeout=timeout,
                        catalog_cfg=cfg,
                    )
                except Exception as exc:  # pragma: no cover - defensive logging
                    logger.warning(
                        "parent_lookup_rebatch_failed",
                        extra={"count": len(chunk), "error": str(exc)},
                    )
                    continue
                result.update(chunk_result)
            return result

        retry_chunk_size = 0
        if allow_rebatch and len(pending) >= _PARENT_LOOKUP_RETRY_THRESHOLD:
            base_chunk_size = max(1, cfg.page_size)
            retry_chunk_size = base_chunk_size // 2
            if retry_chunk_size < 2:
                retry_chunk_size = 0
            elif retry_chunk_size >= base_chunk_size and base_chunk_size > 1:
                retry_chunk_size = base_chunk_size - 1
            retry_chunk_size = min(len(pending), retry_chunk_size)

        remaining: list[str] = []
        if retry_chunk_size and len(pending) > 1:
            for chunk in _chunked(pending, retry_chunk_size):
                batch_attempts += 1
                try:
                    chunk_result = _fetch_parent_catalog_chunk(
                        chunk,
                        client=client,
                        api_cfg=api_cfg,
                        timeout=timeout,
                        catalog_cfg=cfg,
                    )
                except Exception as exc:  # pragma: no cover - defensive logging
                    logger.warning(
                        "parent_lookup_rebatch_failed",
                        extra={"count": len(chunk), "error": str(exc)},
                    )
                    remaining.extend(chunk)
                    continue
                result.update(chunk_result)
                missing = [
                    chembl_id for chembl_id in chunk if chembl_id not in chunk_result
                ]
                if missing:
                    remaining.extend(missing)
        else:
            remaining = [
                chembl_id for chembl_id in pending if chembl_id not in existing
            ]

        outstanding_pool = {
            chembl_id for chembl_id in remaining if chembl_id not in existing
        }
        outstanding: list[str] = []
        if outstanding_pool:
            seen_outstanding: set[str] = set()
            for chembl_id in pending:
                if chembl_id not in outstanding_pool:
                    continue
                if chembl_id in result or chembl_id in seen_outstanding:
                    continue
                seen_outstanding.add(chembl_id)
                outstanding.append(chembl_id)
        if outstanding and not _PARENT_CATALOG_LOCK.locked():
            try:
                load_parent_catalog(
                    client=client,
                    api_cfg=api_cfg,
                    catalog_cfg=cfg,
                    timeout=timeout,
                )
            except Exception as exc:  # pragma: no cover - defensive logging
                logger.warning(
                    "parent_lookup_cache_refresh_failed",
                    extra={"error": str(exc)},
                )
        if outstanding:
            cache_hits = query_parent_catalog(outstanding, cfg)
            if cache_hits:
                result.update(cache_hits)
                outstanding = [
                    chembl_id
                    for chembl_id in outstanding
                    if chembl_id not in cache_hits
                ]

        if outstanding:
            if single_limit is not None and len(outstanding) > single_limit:
                logger.warning(
                    "parent_lookup_single_limit_reached",
                    extra={
                        "limit": single_limit,
                        "remaining": len(outstanding) - single_limit,
                    },
                )
                outstanding = outstanding[:single_limit]

        if outstanding:
            single_requests = len(outstanding)
            max_workers = max(
                1,
                min(
                    len(outstanding),
                    _PARENT_LOOKUP_SINGLE_CONCURRENCY,
                    max(1, api_cfg.rps),
                ),
            )

            resolution_cache: dict[str, tuple[str, str] | None] = {}
            resolution_lock = threading.Lock()
            resolved: dict[str, tuple[str, str] | None] = {}

            def _fetch_with_retry(chembl_id: str) -> tuple[str, str] | None:
                with resolution_lock:
                    if chembl_id in resolution_cache:
                        return resolution_cache[chembl_id]

                pair: tuple[str, str] | None = None
                for attempt in range(1, attempts + 1):
                    try:
                        pair = fetch_parent_for_id(
                            chembl_id,
                            client=client,
                            api_cfg=api_cfg,
                            timeout=timeout,
                            catalog_cfg=cfg,
                        )
                    except Exception as exc:  # pragma: no cover - defensive logging
                        if attempt >= attempts:
                            logger.warning(
                                "parent_lookup_single_failed",
                                extra={"chembl_id": chembl_id, "error": str(exc)},
                            )
                            pair = None
                            break
                        if delay:
                            sleep(delay * (2 ** (attempt - 1)))
                        continue
                    break

                with resolution_lock:
                    resolution_cache[chembl_id] = pair
                return pair

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_map = {
                    executor.submit(_fetch_with_retry, chembl_id): chembl_id
                    for chembl_id in outstanding
                }
                for future in as_completed(future_map):
                    chembl_id = future_map[future]
                    try:
                        resolved[chembl_id] = future.result()
                    except Exception as exc:  # pragma: no cover - defensive logging
                        logger.warning(
                            "parent_lookup_single_failed",
                            extra={"chembl_id": chembl_id, "error": str(exc)},
                        )
                        with resolution_lock:
                            resolution_cache.setdefault(chembl_id, None)
                        resolved[chembl_id] = None

            for chembl_id in outstanding:
                pair = resolved.get(chembl_id)
                if not pair:
                    continue
                child_id, parent_id = pair
                result[child_id] = parent_id
    finally:
        elapsed = perf_counter() - fallback_start
        logger.info(
            "parent_lookup_fallback_summary",
            extra={
                "attempted": len(pending),
                "resolved": len(result),
                "single_ids": single_requests,
                "batch_attempts": batch_attempts,
                "elapsed_seconds": elapsed,
            },
        )

    return result


def fetch_parent_catalog_for(
    ids: Iterable[str],
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    timeout: float | None = None,
    catalog_cfg: MoleculeCatalogCfg | None = None,
) -> dict[str, str]:
    """Fetch parent mappings for *ids* from the ChEMBL API."""

    normalised_ids: list[str] = []
    for value in ids:
        try:
            normalised_ids.append(_normalise_chembl_id(str(value)))
        except ValueError:
            continue

    unique_ids: list[str] = []
    seen: set[str] = set()
    for item in normalised_ids:
        if item in seen:
            continue
        seen.add(item)
        unique_ids.append(item)

    if not unique_ids:
        return {}

    cfg = catalog_cfg or _DEFAULT_CATALOG_CFG
    effective_timeout = timeout if timeout is not None else api_cfg.timeout_read
    result: dict[str, str] = {}
    parentless_filtered = _filters_exclude_parentless(cfg)

    if len(unique_ids) <= _PARENT_LOOKUP_FALLBACK_THRESHOLD:
        missing: list[str] = unique_ids
        use_bulk_lookup = parentless_filtered and len(unique_ids) == 1
        if use_bulk_lookup:
            chunk_result = _fetch_parent_catalog_chunk(
                unique_ids,
                client=client,
                api_cfg=api_cfg,
                timeout=effective_timeout,
                catalog_cfg=cfg,
            )
            result.update(chunk_result)
            missing = [
                chembl_id for chembl_id in unique_ids if chembl_id not in chunk_result
            ]
            if not missing:
                return result
        fallback_result = _fetch_parent_catalog_via_helper(
            missing,
            client=client,
            api_cfg=api_cfg,
            timeout=effective_timeout,
            existing=result,
            catalog_cfg=cfg,
            allow_rebatch=False,
            allow_single_lookup=not parentless_filtered,
        )
        result.update(fallback_result)
        return result

    fallback_candidates: list[str] = []

    chunk_size = max(1, cfg.page_size)
    for chunk in _chunked(unique_ids, chunk_size):
        try:
            chunk_result = _fetch_parent_catalog_chunk(
                chunk,
                client=client,
                api_cfg=api_cfg,
                timeout=effective_timeout,
                catalog_cfg=cfg,
            )
        except Exception as exc:  # pragma: no cover - defensive logging
            logger.warning(
                "parent_lookup_bulk_failed",
                extra={"count": len(chunk), "error": str(exc)},
            )
            fallback_candidates.extend(chunk)
            continue
        result.update(chunk_result)
        if len(chunk_result) != len(chunk):
            missing = [
                chembl_id for chembl_id in chunk if chembl_id not in chunk_result
            ]
            if missing:
                fallback_candidates.extend(missing)

    if fallback_candidates:
        ordered = [
            item for item in dict.fromkeys(fallback_candidates) if item not in result
        ]
        if ordered:
            fallback_result = _fetch_parent_catalog_via_helper(
                ordered,
                client=client,
                api_cfg=api_cfg,
                timeout=effective_timeout,
                existing=result,
                catalog_cfg=cfg,
                allow_single_lookup=not parentless_filtered,
            )
            result.update(fallback_result)

    return result


def fetch_parent_catalog(
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    catalog_cfg: MoleculeCatalogCfg,
    timeout: float | None = None,
) -> dict[str, str]:
    """Fetch the molecule parent catalogue from ChEMBL."""

    url = _build_initial_url(api_cfg, catalog_cfg)
    next_url: str | None = url
    result: dict[str, str] = {}
    effective_timeout = timeout if timeout is not None else api_cfg.timeout_read

    while next_url:
        logger.debug("parent_catalog_page", url=next_url, collected=len(result))
        data = client.request_json(next_url, cfg=api_cfg, timeout=effective_timeout)
        items = (
            data.get("molecules")
            or data.get("molecule")
            or data.get("molecule_parents")
            or data.get("molecule_parent")
            or []
        )
        if not isinstance(items, list):
            logger.warning("unexpected_response_items", extra={"url": next_url})
            items = []
        for item in items:
            if not isinstance(item, Mapping):
                continue
            child = item.get(catalog_cfg.child_field)
            parent = item.get(catalog_cfg.parent_field)
            if not child or not parent:
                continue
            try:
                child_id = _normalise_chembl_id(str(child))
                parent_id = _normalise_chembl_id(str(parent))
            except ValueError:
                continue
            result[child_id] = parent_id
        page_meta = data.get("page_meta") or {}
        next_token = page_meta.get("next")
        next_url = urljoin(api_cfg.chembl_base, next_token) if next_token else None

    return result


def _read_cache(path: Path, catalog_cfg: MoleculeCatalogCfg) -> dict[str, str]:
    if not path.is_file():
        return {}
    suffix = path.suffix.lower()
    if suffix == ".csv":
        try:
            with path.open("r", encoding="utf-8", newline="") as fh:
                reader = csv.DictReader(fh)
                required = {catalog_cfg.child_field, catalog_cfg.parent_field}
                headers = set(reader.fieldnames or ())
                missing = required - headers
                if missing:
                    raise ValueError(
                        "missing columns in parent catalog: "
                        + ", ".join(sorted(missing))
                    )
                parsed: dict[str, str] = {}
                for row in reader:
                    child = row.get(catalog_cfg.child_field)
                    parent = row.get(catalog_cfg.parent_field)
                    if not child or not parent:
                        continue
                    try:
                        child_id = _normalise_chembl_id(str(child))
                        parent_id = _normalise_chembl_id(str(parent))
                    except ValueError:
                        continue
                    parsed[child_id] = parent_id
                return parsed
        except csv.Error as exc:  # pragma: no cover - defensive
            raise ValueError(f"invalid CSV catalog: {path}: {exc}") from exc
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        logger.warning("invalid_catalog_cache", extra={"path": str(path)})
        return {}
    if not isinstance(raw, Mapping):
        logger.warning("unexpected_catalog_cache", extra={"path": str(path)})
        return {}
    result: dict[str, str] = {}
    for child, parent in raw.items():
        if not isinstance(child, str) or not isinstance(parent, str):
            continue
        try:
            child_id = _normalise_chembl_id(child)
            parent_id = _normalise_chembl_id(parent)
        except ValueError:
            continue
        result[child_id] = parent_id
    return result


def _ensure_sqlite_schema(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS parent_catalog (
            child TEXT PRIMARY KEY,
            parent TEXT NOT NULL
        )
        """
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_parent_catalog_parent ON parent_catalog(parent)"
    )


def _write_cache_from_items(
    items: Iterable[tuple[str, str]], catalog_cfg: MoleculeCatalogCfg
) -> None:
    cache_path = catalog_cfg.cache_path
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    entries = list(items)

    if cache_path.suffix.lower() == ".csv":
        with open_atomic(cache_path, encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[catalog_cfg.child_field, catalog_cfg.parent_field],
            )
            writer.writeheader()
            for child, parent in entries:
                writer.writerow(
                    {
                        catalog_cfg.child_field: child,
                        catalog_cfg.parent_field: parent,
                    }
                )
        return

    with open_atomic(cache_path, encoding="utf-8") as fh:
        json.dump(dict(entries), fh, indent=2, sort_keys=True)


def _write_parent_catalog_sqlite(
    items: Iterable[tuple[str, str]],
    catalog_cfg: MoleculeCatalogCfg,
    *,
    replace: bool,
) -> None:
    sqlite_path = catalog_cfg.sqlite_path
    sqlite_path.parent.mkdir(parents=True, exist_ok=True)
    entries = list(items)

    with sqlite3.connect(sqlite_path) as conn:
        _ensure_sqlite_schema(conn)
        if replace:
            conn.execute("DELETE FROM parent_catalog")
        if entries:
            conn.executemany(
                "INSERT OR REPLACE INTO parent_catalog(child, parent) VALUES (?, ?)",
                entries,
            )
        conn.commit()


def _read_sqlite_cache(path: Path) -> tuple[dict[str, str], bool]:
    if not path.is_file():
        return {}, False
    try:
        with sqlite3.connect(path) as conn:
            _ensure_sqlite_schema(conn)
            rows = conn.execute("SELECT child, parent FROM parent_catalog")
            result = {str(child): str(parent) for child, parent in rows}
        return result, True
    except sqlite3.DatabaseError as exc:
        logger.warning(
            "invalid_catalog_sqlite",
            extra={"path": str(path), "error": str(exc)},
        )
        return {}, False


def _ensure_sqlite_cache(catalog_cfg: MoleculeCatalogCfg) -> bool:
    sqlite_path = catalog_cfg.sqlite_path
    if sqlite_path.is_file():
        return True
    cached = _read_cache(catalog_cfg.cache_path, catalog_cfg)
    if not cached:
        return False
    items = sorted((str(child), str(parent)) for child, parent in cached.items())
    _write_parent_catalog_sqlite(items, catalog_cfg, replace=True)
    return sqlite_path.is_file()


def query_parent_catalog(
    children: Iterable[str], catalog_cfg: MoleculeCatalogCfg
) -> dict[str, str]:
    """Return parent mappings for *children* from the SQLite cache."""

    unique: list[str] = []
    seen: set[str] = set()
    for value in children:
        try:
            child_id = _normalise_chembl_id(str(value))
        except ValueError:
            continue
        if child_id in seen:
            continue
        seen.add(child_id)
        unique.append(child_id)

    if not unique:
        return {}

    if not _ensure_sqlite_cache(catalog_cfg):
        return {}

    sqlite_path = catalog_cfg.sqlite_path
    result: dict[str, str] = {}

    try:
        with sqlite3.connect(sqlite_path) as conn:
            _ensure_sqlite_schema(conn)
            for chunk in _chunked(unique, _SQLITE_VARIABLE_LIMIT):
                if not chunk:
                    continue
                placeholders = ",".join("?" for _ in chunk)
                cursor = conn.execute(
                    f"SELECT child, parent FROM parent_catalog WHERE child IN ({placeholders})",
                    chunk,
                )
                result.update({str(child): str(parent) for child, parent in cursor})
    except sqlite3.DatabaseError as exc:
        logger.warning(
            "invalid_catalog_sqlite",
            extra={"path": str(sqlite_path), "error": str(exc)},
        )
        return {}

    return result


def write_parent_catalog_cache(
    catalog: Mapping[str, str], catalog_cfg: MoleculeCatalogCfg
) -> None:
    """Persist *catalog* to disk using the configured cache format."""

    sorted_items = sorted(
        (str(child), str(parent)) for child, parent in catalog.items()
    )
    _write_cache_from_items(sorted_items, catalog_cfg)
    _write_parent_catalog_sqlite(sorted_items, catalog_cfg, replace=True)


def update_parent_catalog_cache(
    catalog: Mapping[str, str], catalog_cfg: MoleculeCatalogCfg
) -> None:
    """Update the existing cache with *catalog* entries."""

    if not catalog:
        return

    sorted_items = sorted(
        (str(child), str(parent)) for child, parent in catalog.items()
    )
    _write_parent_catalog_sqlite(sorted_items, catalog_cfg, replace=False)
    data, ok = _read_sqlite_cache(catalog_cfg.sqlite_path)
    if ok:
        _write_cache_from_items(sorted(data.items()), catalog_cfg)


def load_parent_catalog(
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    catalog_cfg: MoleculeCatalogCfg,
    timeout: float | None = None,
    force_refresh: bool = False,
) -> dict[str, str]:
    """Return the molecule parent catalogue, using the on-disk cache if present."""

    with _PARENT_CATALOG_LOCK:
        cache_path = catalog_cfg.cache_path
        sqlite_path = catalog_cfg.sqlite_path
        if not force_refresh:
            sqlite_data, sqlite_ok = _read_sqlite_cache(sqlite_path)
            if sqlite_ok:
                return sqlite_data
            cached = _read_cache(cache_path, catalog_cfg)
            if cached:
                items = sorted(cached.items())
                _write_parent_catalog_sqlite(items, catalog_cfg, replace=True)
                return cached

        result = fetch_parent_catalog(
            client=client, api_cfg=api_cfg, catalog_cfg=catalog_cfg, timeout=timeout
        )
        write_parent_catalog_cache(result, catalog_cfg)
        return result
