"""Utilities for retrieving and caching the ChEMBL molecule parent catalogue."""

from __future__ import annotations

import csv
import json
import sqlite3
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import Mapping
from urllib.parse import urlencode, urljoin

from .chembl_client import ChemblClient
from .config import ApiCfg, MoleculeCatalogCfg
from .log import logger

__all__ = [
    "ParentCatalogStore",
    "fetch_parent_catalog",
    "load_parent_catalog",
]

_SQLITE_PARAMETER_LIMIT = 999


class ParentCatalogStore:
    """Read-only view over the cached parent catalogue stored in SQLite."""

    def __init__(self, path: Path, *, refreshed: bool) -> None:
        self._path = path
        self._refreshed = refreshed

    @property
    def path(self) -> Path:
        """Return the backing SQLite database path."""

        return self._path

    @property
    def was_refreshed(self) -> bool:
        """Whether the catalogue was refreshed from the remote endpoint."""

        return self._refreshed

    def __bool__(self) -> bool:
        return _sqlite_has_rows(self._path)

    def lookup(self, children: Sequence[str]) -> dict[str, str]:
        """Return a mapping for ``children`` found in the cache."""

        if not children:
            return {}

        # ``children`` are expected to be normalised already, yet we only
        # retain non-empty identifiers to avoid unnecessary queries.
        filtered = list(dict.fromkeys(child for child in children if child))
        if not filtered:
            return {}

        result: dict[str, str] = {}
        with sqlite3.connect(self._path) as connection:
            for chunk in _chunked(filtered, _SQLITE_PARAMETER_LIMIT):
                placeholders = ",".join("?" for _ in chunk)
                query = (
                    f"SELECT child, parent FROM parent_catalog "
                    f"WHERE child IN ({placeholders})"
                )
                for child, parent in connection.execute(query, tuple(chunk)):
                    result[str(child)] = str(parent)
        return result


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


def _write_sqlite_cache(path: Path, data: Mapping[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(path) as connection:
        connection.execute("DROP TABLE IF EXISTS parent_catalog")
        connection.execute(
            "CREATE TABLE parent_catalog (child TEXT PRIMARY KEY, parent TEXT NOT NULL)"
        )
        rows = sorted((str(child), str(parent)) for child, parent in data.items())
        connection.executemany(
            "INSERT INTO parent_catalog (child, parent) VALUES (?, ?)",
            rows,
        )
        connection.execute(
            "CREATE INDEX IF NOT EXISTS parent_catalog_parent_idx"
            " ON parent_catalog (parent)"
        )


def _sqlite_has_rows(path: Path) -> bool:
    if not path.is_file():
        return False
    with sqlite3.connect(path) as connection:
        cursor = connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='parent_catalog'"
        )
        if cursor.fetchone() is None:
            return False
        cursor = connection.execute("SELECT 1 FROM parent_catalog LIMIT 1")
        return cursor.fetchone() is not None


def _chunked(values: Iterable[str], limit: int) -> Iterator[list[str]]:
    if limit <= 0:
        raise ValueError("chunk size must be positive")
    chunk: list[str] = []
    for value in values:
        chunk.append(value)
        if len(chunk) >= limit:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def load_parent_catalog(
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    catalog_cfg: MoleculeCatalogCfg,
    timeout: float | None = None,
    force_refresh: bool = False,
) -> ParentCatalogStore:
    """Ensure the parent catalogue cache exists and return a SQLite-backed view."""

    sqlite_path = catalog_cfg.sqlite_path
    refreshed = False

    if not force_refresh and _sqlite_has_rows(sqlite_path):
        return ParentCatalogStore(sqlite_path, refreshed=refreshed)

    if not force_refresh:
        cached = _read_cache(catalog_cfg.cache_path, catalog_cfg)
        if cached:
            _write_sqlite_cache(sqlite_path, cached)
            return ParentCatalogStore(sqlite_path, refreshed=refreshed)

    result = fetch_parent_catalog(
        client=client, api_cfg=api_cfg, catalog_cfg=catalog_cfg, timeout=timeout
    )
    _write_sqlite_cache(sqlite_path, result)
    refreshed = True
    return ParentCatalogStore(sqlite_path, refreshed=refreshed)
