"""Utilities for retrieving and caching the ChEMBL molecule parent catalogue."""

from __future__ import annotations

import csv
import json
import sqlite3
from collections.abc import Iterable
from pathlib import Path
from typing import Mapping
from urllib.parse import urlencode, urljoin

from .chembl_client import ChemblClient, _chunked
from .config import ApiCfg, MoleculeCatalogCfg
from .log import logger

_DEFAULT_CATALOG_CFG = MoleculeCatalogCfg()
_PARENT_LOOKUP_FIELDS = tuple(_DEFAULT_CATALOG_CFG.fields or ()) or (
    _DEFAULT_CATALOG_CFG.child_field,
    _DEFAULT_CATALOG_CFG.parent_field,
)
_PARENT_LOOKUP_CHILD_FIELD = _DEFAULT_CATALOG_CFG.child_field
_PARENT_LOOKUP_PARENT_FIELD = _DEFAULT_CATALOG_CFG.parent_field
_PARENT_LOOKUP_CHUNK_SIZE = _DEFAULT_CATALOG_CFG.page_size
_SQLITE_VARIABLE_LIMIT = 900

__all__ = [
    "fetch_parent_catalog",
    "fetch_parent_catalog_for",
    "load_parent_catalog",
    "query_parent_catalog",
    "update_parent_catalog_cache",
    "write_parent_catalog_cache",
]


def _normalise_chembl_id(value: str) -> str:
    """Return *value* normalised as an upper-case ChEMBL identifier."""

    normalised = value.strip().upper()
    if not normalised:
        raise ValueError("empty identifier")
    return normalised


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


def fetch_parent_catalog_for(
    ids: Iterable[str],
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    timeout: float | None = None,
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

    base = api_cfg.chembl_base.rstrip("/")
    effective_timeout = timeout if timeout is not None else api_cfg.timeout_read
    result: dict[str, str] = {}

    for chunk in _chunked(unique_ids, _PARENT_LOOKUP_CHUNK_SIZE):
        params = {
            "format": "json",
            "limit": str(len(chunk)),
            f"{_PARENT_LOOKUP_CHILD_FIELD}__in": ",".join(chunk),
            "fields": ",".join(_PARENT_LOOKUP_FIELDS),
        }
        for key, value in _DEFAULT_CATALOG_CFG.filters.items():
            params[key] = value
        url = f"{base}/molecule.json?{urlencode(params)}"
        data = client.request_json(url, cfg=api_cfg, timeout=effective_timeout)
        items = (
            data.get("molecules")
            or data.get("molecule")
            or data.get("molecule_parents")
            or data.get("molecule_parent")
            or []
        )
        if not isinstance(items, list):
            logger.warning("unexpected_response_items", extra={"url": url})
            continue
        allowed = set(chunk)
        for item in items:
            if not isinstance(item, Mapping):
                continue
            child = item.get(_PARENT_LOOKUP_CHILD_FIELD)
            parent = item.get(_PARENT_LOOKUP_PARENT_FIELD)
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
                result: dict[str, str] = {}
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
                    result[child_id] = parent_id
                return result
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
        with cache_path.open("w", encoding="utf-8", newline="") as fh:
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

    cache_path.write_text(
        json.dumps(dict(entries), indent=2, sort_keys=True),
        encoding="utf-8",
    )


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
