"""Low-level PubChem client utilities."""

from __future__ import annotations

from dataclasses import dataclass
from threading import Lock
from time import monotonic
from typing import Any, cast
from urllib.parse import quote

import requests
from cachetools import TTLCache
from requests import Session

from ..config import ApiCfg, PubChemCfg, RetryCfg, session_with_retry
from ..log import logger
from ..rate_limiter import get_limiter, sleep

__all__ = [
    "Properties",
    "init_session",
    "make_request",
    "url_encode",
    "get_cid_from_smiles",
    "get_cid_from_inchi",
    "get_cid_from_inchikey",
    "get_cid",
    "get_all_cid",
    "get_standard_name",
    "get_properties",
    "validate_cid",
]


# Cache is initialised lazily to allow configuration of the TTL via
# :class:`PubChemCfg`. The cache is recreated when the TTL changes.
_CACHE: TTLCache[str, dict[str, Any]] | None = None
_CACHE_LOCK = Lock()

# Shared session with placeholder user agent; production code should call
# :func:`init_session` to supply real contact details.
_session: Session = session_with_retry(
    ApiCfg(user_agent="chembl-da/0.1 (mailto:contact@example.org)"), RetryCfg()
)


def init_session(api: ApiCfg, retry: RetryCfg) -> None:
    """Initialise the shared HTTP session."""

    global _session
    _session = session_with_retry(api, retry)


def url_encode(text: str) -> str:
    """Return ``text`` URL-encoded for safe usage in HTTP requests."""

    return quote(text, safe="")


def _cids_from_identifier_list(data: dict[str, Any]) -> list[str]:
    """Extract CIDs from a JSON ``IdentifierList`` structure."""

    return [str(cid) for cid in data.get("IdentifierList", {}).get("CID", [])]


def get_cid_from_smiles(smiles: str, cfg: PubChemCfg) -> str | None:
    """Retrieve PubChem CID(s) for a SMILES string."""

    safe_smiles = url_encode(smiles)
    base = cfg.base.rstrip("/")
    url = f"{base}/compound/smiles/{safe_smiles}/cids/JSON"
    response = make_request(url, cfg)
    if not response:
        return None
    cids = _cids_from_identifier_list(response)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_cid_from_inchi(inchi: str, cfg: PubChemCfg) -> str | None:
    """Retrieve PubChem CID(s) for an InChI string."""

    safe_inchi = url_encode(inchi)
    base = cfg.base.rstrip("/")
    url = f"{base}/compound/inchi/{safe_inchi}/cids/JSON"
    response = make_request(url, cfg)
    if not response:
        return None
    cids = _cids_from_identifier_list(response)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_cid_from_inchikey(inchikey: str, cfg: PubChemCfg) -> str | None:
    """Retrieve PubChem CID(s) for an InChIKey."""

    safe_inchikey = url_encode(inchikey)
    base = cfg.base.rstrip("/")
    url = f"{base}/compound/inchikey/{safe_inchikey}/cids/JSON"
    response = make_request(url, cfg)
    if not response:
        return None
    cids = _cids_from_identifier_list(response)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def make_request(url: str, cfg: PubChemCfg) -> dict[str, Any] | None:
    """Make an HTTP GET request and return parsed JSON."""

    global _CACHE
    with _CACHE_LOCK:
        if _CACHE is None or _CACHE.ttl != cfg.cache_ttl:
            _CACHE = TTLCache(maxsize=1024, ttl=cfg.cache_ttl)
        cache = _CACHE
        cached = cache.get(url) if cache is not None else None
    if cached is not None:
        logger.info("cache_hit", url=url, rps=cfg.rps, status="hit")
        return cast(dict[str, Any], cached)
    logger.info("cache_miss", url=url, rps=cfg.rps, status="miss")

    attempts = max(cfg.retries, 1)
    backoff_delay = cfg.backoff_initial_seconds
    deadline: float | None = None
    if cfg.timeout_seconds > 0:
        deadline = monotonic() + cfg.timeout_seconds

    for attempt in range(1, attempts + 1):
        if deadline is not None and monotonic() >= deadline:
            logger.warning("request_timeout", url=url, attempt=attempt, rps=cfg.rps)
            logger.info("request_fail", url=url, status=None, rps=cfg.rps)
            return None
        event = "request_start" if attempt == 1 else "request_retry"
        logger.info(event, url=url, attempt=attempt, rps=cfg.rps)
        get_limiter("pubchem", cfg.rps, cfg.burst).acquire()
        try:
            response = _session.get(
                url, timeout=(cfg.timeout_connect, cfg.timeout_read)
            )
        except requests.RequestException as exc:  # pragma: no cover - network
            if attempt >= attempts:
                logger.error(
                    "request_error",
                    url=url,
                    error=str(exc),
                    attempt=attempt,
                    rps=cfg.rps,
                )
                logger.info(
                    "request_fail",
                    url=url,
                    status=None,
                    rps=cfg.rps,
                )
                return None
            if cfg.delay > 0:
                sleep(cfg.delay)
            continue

        status = response.status_code
        if status == 404:
            logger.info("request_not_found", url=url, status=status, rps=cfg.rps)
            logger.info("request_fail", url=url, status=status, rps=cfg.rps)
            return None
        if status == 429 or 500 <= status < 600:
            event_name = (
                "request_rate_limited" if status == 429 else "request_server_error"
            )
            logger.warning(
                event_name,
                url=url,
                status=status,
                attempt=attempt,
                rps=cfg.rps,
            )
            if attempt >= attempts:
                logger.info("request_fail", url=url, status=status, rps=cfg.rps)
                return None
            delay = backoff_delay if backoff_delay > 0 else cfg.delay
            if delay > 0:
                sleep(delay)
                backoff_delay = delay * 2
            continue
        if status >= 400:
            logger.warning(
                "request_unexpected_status",
                url=url,
                status=status,
                rps=cfg.rps,
            )
            logger.info("request_fail", url=url, status=status, rps=cfg.rps)
            return None

        try:
            response.raise_for_status()
            data = cast(dict[str, Any], response.json())
        except requests.RequestException as exc:  # pragma: no cover - network
            if attempt >= attempts:
                logger.error(
                    "request_error",
                    url=url,
                    error=str(exc),
                    attempt=attempt,
                    rps=cfg.rps,
                )
                logger.info("request_fail", url=url, status=status, rps=cfg.rps)
                return None
            if cfg.delay > 0:
                sleep(cfg.delay)
            continue
        except ValueError:
            logger.warning(
                "response_not_json",
                url=url,
                status=status,
                rps=cfg.rps,
            )
            logger.info("request_fail", url=url, status=status, rps=cfg.rps)
            return None

        logger.info(
            "request_ok",
            url=url,
            status=status,
            rps=cfg.rps,
        )
        with _CACHE_LOCK:
            cache = _CACHE
            if cache is None or cache.ttl != cfg.cache_ttl:
                cache = _CACHE = TTLCache(maxsize=1024, ttl=cfg.cache_ttl)
            cache[url] = data
        logger.info("cache_set", url=url, rps=cfg.rps)
        return data

    return None


def _extract_cids(bindings: list[dict[str, Any]]) -> list[str]:
    """Extract CIDs from API bindings."""

    cids: list[str] = []
    for item in bindings:
        cid_field = item.get("cid")
        if isinstance(cid_field, dict):
            cid_value = cid_field.get("value", "")
        else:
            cid_value = str(cid_field)
        cid_value = cid_value.replace(
            "http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID", ""
        )
        if cid_value:
            cids.append(cid_value)
    return cids


def get_cid(compound_name: str, cfg: PubChemCfg) -> str | None:
    """Retrieve PubChem CID(s) for *compound_name* (exact match)."""

    safe_name = url_encode(compound_name)
    rdf_base = cfg.base.rstrip("/").rsplit("/", 1)[0] + "/rdf"
    url = f"{rdf_base}/query?graph=synonym&return=cid&format=json&name={safe_name}"
    response = make_request(url, cfg)
    if not response:
        return None
    bindings = response.get("results", {}).get("bindings", [])
    cids = _extract_cids(bindings)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_all_cid(compound_name: str, cfg: PubChemCfg) -> str | None:
    """Retrieve PubChem CID(s) for *compound_name* (partial match)."""

    safe_name = url_encode(compound_name)
    rdf_base = cfg.base.rstrip("/").rsplit("/", 1)[0] + "/rdf"
    url = f"{rdf_base}/query?graph=synonym&return=cid&format=json&name={safe_name}&contain=true"
    response = make_request(url, cfg)
    if not response:
        return None
    bindings = response.get("results", {}).get("bindings", [])
    cids = _extract_cids(bindings)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def validate_cid(cid: str) -> str | None:
    """Validate PubChem CID."""

    if cid in {"", "0", "-1"}:
        return None
    return cid


def get_standard_name(cid: str, cfg: PubChemCfg) -> str | None:
    """Retrieve the standard compound name for a given CID."""

    validated = validate_cid(cid)
    if not validated:
        return None
    base = cfg.base.rstrip("/")
    url = f"{base}/compound/cid/{validated}/description/JSON"
    response = make_request(url, cfg)
    if not response:
        return None
    info = response.get("InformationList", {}).get("Information", [])
    if not info:
        return None
    return cast(str | None, info[0].get("Title"))


@dataclass
class Properties:
    """Chemical properties for a PubChem compound."""

    IUPACName: str | None
    MolecularFormula: str | None
    iSMILES: str | None
    cSMILES: str | None
    InChI: str | None
    InChIKey: str | None


def get_properties(cid: str, cfg: PubChemCfg) -> Properties:
    """Retrieve chemical properties for a compound by CID."""

    validated = validate_cid(cid)
    if not validated:
        return Properties(None, None, None, None, None, None)
    base = cfg.base.rstrip("/")
    url = (
        f"{base}/compound/cid/{validated}/property/MolecularFormula,IUPACName,IsomericSMILES,"
        "CanonicalSMILES,InChI,InChIKey/JSON"
    )
    response = make_request(url, cfg)
    if not response:
        return Properties(None, None, None, None, None, None)
    props = response.get("PropertyTable", {}).get("Properties", [])
    if not props:
        return Properties(None, None, None, None, None, None)
    item = props[0]
    return Properties(
        cast(str | None, item.get("IUPACName")),
        cast(str | None, item.get("MolecularFormula")),
        cast(str | None, item.get("IsomericSMILES")),
        cast(str | None, item.get("CanonicalSMILES")),
        cast(str | None, item.get("InChI")),
        cast(str | None, item.get("InChIKey")),
    )
