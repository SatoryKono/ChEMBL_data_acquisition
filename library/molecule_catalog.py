"""Utilities for retrieving and caching the ChEMBL molecule parent catalogue."""

from __future__ import annotations

import csv
import json
import sqlite3
from collections.abc import Iterable, Mapping as ABCMapping
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
_SQLITE_TABLE = "parent_catalog"
_SQLITE_CREATE_TABLE = (
    f"CREATE TABLE IF NOT EXISTS {_SQLITE_TABLE} ("
    "child TEXT PRIMARY KEY,"
    "parent TEXT NOT NULL"
    ")"
)
_SQLITE_PARENT_INDEX = (
    f"CREATE INDEX IF NOT EXISTS ix_{_SQLITE_TABLE}_parent ON {_SQLITE_TABLE} (parent)"
)
_SQLITE_MAX_PARAMETERS = 900


class ParentCatalog:
    """SQLite-backed view of the parent molecule catalogue."""

    def __init__(
        self,
        sqlite_path: Path,
        *,
        has_rows: bool,
        refreshed: bool,
    ) -> None:
        self._sqlite_path = sqlite_path
        self._has_rows = has_rows
        self.refreshed = refreshed

    @property
    def sqlite_path(self) -> Path:
        """Return the backing SQLite database path."""

        return self._sqlite_path

    def __bool__(self) -> bool:  # pragma: no cover - trivial
        return self._has_rows

    def lookup_many(self, children: Iterable[str]) -> dict[str, str]:
        """Return parent mappings for ``children`` from the SQLite cache."""

        normalised: list[str] = []
        seen: set[str] = set()
        for value in children:
            try:
                child_id = _normalise_chembl_id(str(value))
            except ValueError:
                continue
            if child_id in seen:
                continue
            seen.add(child_id)
            normalised.append(child_id)

        if not normalised:
            return {}

        result: dict[str, str] = {}
        with sqlite3.connect(self._sqlite_path) as conn:
            conn.row_factory = sqlite3.Row
            for chunk in _chunked(normalised, _SQLITE_MAX_PARAMETERS):
                placeholders = ",".join("?" for _ in chunk)
                query = (
                    f"SELECT child, parent FROM {_SQLITE_TABLE} "
                    f"WHERE child IN ({placeholders})"
                )
                cursor = conn.execute(query, chunk)
                result.update({row["child"]: row["parent"] for row in cursor})
        return result

    def upsert_many(self, mapping: Mapping[str, str]) -> None:
        """Insert or update ``mapping`` in the SQLite cache."""

        rows: list[tuple[str, str]] = []
        for child, parent in mapping.items():
            try:
                child_id = _normalise_chembl_id(str(child))
                parent_id = _normalise_chembl_id(str(parent))
            except ValueError:
                continue
            rows.append((child_id, parent_id))

        if not rows:
            return

        with sqlite3.connect(self._sqlite_path) as conn:
            conn.execute("BEGIN")
            conn.executemany(
                f"INSERT INTO {_SQLITE_TABLE} (child, parent) VALUES (?, ?) "
                "ON CONFLICT(child) DO UPDATE SET parent=excluded.parent",
                rows,
            )
            conn.commit()
        self._has_rows = True

    def to_dict(self) -> dict[str, str]:
        """Return the entire parent catalogue as a dictionary."""

        with sqlite3.connect(self._sqlite_path) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(f"SELECT child, parent FROM {_SQLITE_TABLE}")
            return {row["child"]: row["parent"] for row in cursor}

__all__ = [
    "ParentCatalog",
    "fetch_parent_catalog",
    "fetch_parent_catalog_for",
    "load_parent_catalog",
]


def _normalise_chembl_id(value: str) -> str:
    """Return *value* normalised as an upper-case ChEMBL identifier."""

    normalised = value.strip().upper()
    if not normalised:
        raise ValueError("empty identifier")
    return normalised


def _ensure_sqlite_schema(sqlite_path: Path) -> None:
    sqlite_path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(sqlite_path) as conn:
        conn.execute(_SQLITE_CREATE_TABLE)
        conn.execute(_SQLITE_PARENT_INDEX)
        conn.commit()


def _sqlite_has_rows(sqlite_path: Path) -> bool:
    with sqlite3.connect(sqlite_path) as conn:
        cursor = conn.execute(f"SELECT 1 FROM {_SQLITE_TABLE} LIMIT 1")
        return cursor.fetchone() is not None


def _normalise_catalog_mapping(mapping: Mapping[str, str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for child, parent in mapping.items():
        try:
            child_id = _normalise_chembl_id(str(child))
            parent_id = _normalise_chembl_id(str(parent))
        except ValueError:
            continue
        result[child_id] = parent_id
    return result


def _replace_catalog(sqlite_path: Path, mapping: Mapping[str, str]) -> bool:
    normalised = _normalise_catalog_mapping(mapping)
    with sqlite3.connect(sqlite_path) as conn:
        conn.execute("BEGIN")
        conn.execute(f"DELETE FROM {_SQLITE_TABLE}")
        if normalised:
            conn.executemany(
                f"INSERT INTO {_SQLITE_TABLE} (child, parent) VALUES (?, ?)",
                list(normalised.items()),
            )
        conn.commit()
    return bool(normalised)


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
            if not isinstance(item, ABCMapping):
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
            if not isinstance(item, ABCMapping):
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
                return _normalise_catalog_mapping(result)
        except csv.Error as exc:  # pragma: no cover - defensive
            raise ValueError(f"invalid CSV catalog: {path}: {exc}") from exc
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        logger.warning("invalid_catalog_cache", extra={"path": str(path)})
        return {}
    if not isinstance(raw, ABCMapping):
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


def load_parent_catalog(
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    catalog_cfg: MoleculeCatalogCfg,
    timeout: float | None = None,
    force_refresh: bool = False,
) -> ParentCatalog:
    """Ensure the parent catalogue cache exists and return a SQLite handle."""

    sqlite_path = catalog_cfg.sqlite_path
    cache_path = catalog_cfg.cache_path
    _ensure_sqlite_schema(sqlite_path)

    refreshed = False

    if force_refresh:
        data = fetch_parent_catalog(
            client=client,
            api_cfg=api_cfg,
            catalog_cfg=catalog_cfg,
            timeout=timeout,
        )
        has_rows = _replace_catalog(sqlite_path, data)
        refreshed = True
        return ParentCatalog(sqlite_path, has_rows=has_rows, refreshed=refreshed)

    has_rows = _sqlite_has_rows(sqlite_path)
    if not has_rows:
        cached = _read_cache(cache_path, catalog_cfg)
        if cached:
            has_rows = _replace_catalog(sqlite_path, cached)
        else:
            data = fetch_parent_catalog(
                client=client,
                api_cfg=api_cfg,
                catalog_cfg=catalog_cfg,
                timeout=timeout,
            )
            has_rows = _replace_catalog(sqlite_path, data)
            refreshed = True

    return ParentCatalog(sqlite_path, has_rows=has_rows, refreshed=refreshed)
