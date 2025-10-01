"""Low-level PubChem client utilities."""

from __future__ import annotations

import json

from dataclasses import dataclass
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
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
    "get_session",
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


def _config_signature(api: ApiCfg, retry: RetryCfg) -> str:
    """Return a stable signature for the session configuration."""

    api_dump = api.model_dump(mode="json")
    retry_dump = retry.model_dump(mode="json")
    return json.dumps({"api": api_dump, "retry": retry_dump}, sort_keys=True)


@dataclass(frozen=True)
class _CacheEntry:
    """Cache entry capturing both payload and outcome."""

    payload: dict[str, Any] | None
    outcome: str
    details: dict[str, Any] | None = None

    @property
    def is_hit(self) -> bool:
        return self.payload is not None


def _retry_after_seconds(value: str | None, *, now: datetime | None = None) -> float | None:
    """Return the delay in seconds represented by ``Retry-After`` header ``value``."""

    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        seconds = float(text)
    except ValueError:
        try:
            parsed = parsedate_to_datetime(text)
        except (TypeError, ValueError):
            return None
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=timezone.utc)
        reference = now or datetime.now(timezone.utc)
        seconds = (parsed - reference).total_seconds()
    if seconds < 0:
        return 0.0
    return seconds


# Cache is initialised lazily to allow configuration of the TTL via
# :class:`PubChemCfg`. The cache is recreated when the TTL changes.
_CACHE: TTLCache[str, _CacheEntry] | None = None
_CACHE_LOCK = Lock()

_SESSION_LOCK = Lock()
_DEFAULT_API_CFG = ApiCfg(user_agent="chembl-da/1.0 (mailto:chembl-data@ebi.ac.uk)")
_DEFAULT_RETRY_CFG = RetryCfg()
_SESSION_CFG: tuple[ApiCfg, RetryCfg] = (_DEFAULT_API_CFG, _DEFAULT_RETRY_CFG)
_SESSION_SIGNATURE = _config_signature(*_SESSION_CFG)
_SESSION_INITIALISED = False
_session: Session | None = session_with_retry(*_SESSION_CFG)

_PLACEHOLDER_CONTACT = "contact@example.org"
_USER_AGENT_ERROR = (
    "PubChem client requires a custom User-Agent; "
    "call init_session with contact details before making requests."
)


def _has_contact_details(user_agent: str | None) -> bool:
    """Return ``True`` when *user_agent* contains usable contact details."""

    if not user_agent:
        return False
    lowered = user_agent.casefold()
    if _PLACEHOLDER_CONTACT in lowered:
        return False
    return "@" in lowered


def init_session(api: ApiCfg, retry: RetryCfg) -> None:
    """Initialise the shared HTTP session."""

    global _SESSION_CFG, _SESSION_SIGNATURE, _SESSION_INITIALISED, _session
    signature = _config_signature(api, retry)
    old_session: Session | None = None
    with _SESSION_LOCK:
        old_session = _session
        _session = None
        _SESSION_CFG = (api, retry)
        _SESSION_SIGNATURE = signature
        _SESSION_INITIALISED = True
    if old_session is not None:
        old_session.close()


def get_session(cfg: ApiCfg | None = None) -> Session:
    """Return the shared HTTP session, creating or refreshing as needed."""

    global _SESSION_CFG, _SESSION_SIGNATURE, _SESSION_INITIALISED, _session
    old_session: Session | None = None
    with _SESSION_LOCK:
        current_api, current_retry = _SESSION_CFG
        target_api = cfg or current_api
        signature = _config_signature(target_api, current_retry)
        needs_refresh = signature != _SESSION_SIGNATURE or _session is None
        _SESSION_CFG = (target_api, current_retry)
        if needs_refresh:
            if not _SESSION_INITIALISED and not _has_contact_details(
                target_api.user_agent
            ):
                raise ValueError(_USER_AGENT_ERROR)
            new_session = session_with_retry(target_api, current_retry)
            old_session = _session
            _session = new_session
            _SESSION_SIGNATURE = signature
            if not _SESSION_INITIALISED:
                _SESSION_INITIALISED = True
        session = _session
    if old_session is not None:
        old_session.close()
    if session is None:
        raise RuntimeError("Failed to initialise PubChem session")
    if not _SESSION_INITIALISED:
        user_agent = session.headers.get("User-Agent")
        if not _has_contact_details(user_agent):
            raise ValueError(_USER_AGENT_ERROR)
        _SESSION_INITIALISED = True
    return session


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
        cache = _ensure_cache(cfg.cache_ttl)
        cached = cache.get(url) if cache is not None else None
    if cached is not None:
        if cached.is_hit:
            logger.debug("cache_hit", url=url, rps=cfg.rps, status="hit")
            return cast(dict[str, Any], cached.payload)
        miss_details = cached.details or {}
        logger.debug(
            "cache_hit",
            url=url,
            rps=cfg.rps,
            status="miss",
            outcome=cached.outcome,
            **miss_details,
        )
        return None
    logger.debug("cache_miss", url=url, rps=cfg.rps, status="miss")

    total_attempts = cfg.retries + 1
    if total_attempts <= 0:
        total_attempts = 1
    backoff_delay = cfg.backoff_initial_seconds
    deadline: float | None = None
    if cfg.timeout_seconds > 0:
        deadline = monotonic() + cfg.timeout_seconds
    last_failure_details: dict[str, Any] | None = None

    def details_excluding(*keys: str) -> dict[str, Any]:
        if not last_failure_details:
            return {}
        return {
            key: value
            for key, value in last_failure_details.items()
            if key not in keys
        }

    for attempt in range(1, total_attempts + 1):
        if deadline is not None and monotonic() >= deadline:
            logger.warning(
                "request_timeout",
                url=url,
                attempt=attempt,
                total_attempts=total_attempts,
                rps=cfg.rps,
                **(last_failure_details or {}),
            )
            _store_cache_miss(url, cfg, "timeout", last_failure_details)
            logger.debug(
                "request_fail",
                url=url,
                status=None,
                total_attempts=total_attempts,
                rps=cfg.rps,
                **details_excluding("status"),
            )
            return None
        event = "request_start" if attempt == 1 else "request_retry"
        logger.debug(
            event,
            url=url,
            attempt=attempt,
            total_attempts=total_attempts,
            rps=cfg.rps,
        )
        get_limiter("pubchem", cfg.rps, cfg.burst).acquire()
        try:
            session = get_session()
            with session.get(
                url, timeout=(cfg.timeout_connect, cfg.timeout_read)
            ) as response:
                status = response.status_code
                if status == 404:
                    last_failure_details = {"reason": "not_found", "status": status}
                    logger.info(
                        "request_not_found",
                        url=url,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    _store_cache_miss(url, cfg, "not_found", last_failure_details)
                    logger.debug(
                        "request_fail",
                        url=url,
                        total_attempts=total_attempts,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    return None
                if status == 429 or 500 <= status < 600:
                    retry_after_header = response.headers.get("Retry-After")
                    retry_after = _retry_after_seconds(retry_after_header)
                    reason = "rate_limited" if status == 429 else "server_error"
                    last_failure_details = {"reason": reason, "status": status}
                    if retry_after is not None:
                        last_failure_details["retry_after"] = retry_after
                    event_name = (
                        "request_rate_limited"
                        if status == 429
                        else "request_server_error"
                    )
                    logger.warning(
                        event_name,
                        url=url,
                        status=status,
                        attempt=attempt,
                        total_attempts=total_attempts,
                        rps=cfg.rps,
                        **({"retry_after": retry_after}
                          if retry_after is not None
                          else {}),
                    )
                    if attempt >= total_attempts:
                        logger.debug(
                            "request_fail",
                            url=url,
                            total_attempts=total_attempts,
                            rps=cfg.rps,
                            status=status,
                            **details_excluding("status"),
                        )
                        return None
                    if retry_after is not None:
                        delay = retry_after
                        backoff_delay = delay * 2 if delay > 0 else backoff_delay
                    else:
                        delay = backoff_delay if backoff_delay > 0 else cfg.delay
                        if delay > 0:
                            backoff_delay = delay * 2
                    if delay > 0:
                        sleep(delay)
                    continue
                if status >= 400:
                    last_failure_details = {
                        "reason": "unexpected_status",
                        "status": status,
                    }
                    logger.warning(
                        "request_unexpected_status",
                        url=url,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    logger.debug(
                        "request_fail",
                        url=url,
                        total_attempts=total_attempts,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    return None

                try:
                    response.raise_for_status()
                    data = cast(dict[str, Any], response.json())
                except requests.RequestException as exc:  # pragma: no cover - network
                    last_failure_details = {
                        "reason": "response_error",
                        "status": status,
                        "error": str(exc),
                    }
                    if attempt >= total_attempts:
                        logger.error(
                            "request_error",
                            url=url,
                            error=str(exc),
                            attempt=attempt,
                            total_attempts=total_attempts,
                            rps=cfg.rps,
                        )
                        logger.debug(
                            "request_fail",
                            url=url,
                            total_attempts=total_attempts,
                            rps=cfg.rps,
                            status=status,
                            **details_excluding("status"),
                        )
                        return None
                    if cfg.delay > 0:
                        sleep(cfg.delay)
                    continue
                except ValueError:
                    last_failure_details = {
                        "reason": "response_not_json",
                        "status": status,
                    }
                    logger.warning(
                        "response_not_json",
                        url=url,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    logger.debug(
                        "request_fail",
                        url=url,
                        total_attempts=total_attempts,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    return None

                logger.debug(
                    "request_ok",
                    url=url,
                    status=status,
                    rps=cfg.rps,
                )
                with _CACHE_LOCK:
                    cache = _ensure_cache(cfg.cache_ttl)
                    cache[url] = _CacheEntry(payload=data, outcome="hit")
                logger.info("cache_set", url=url, rps=cfg.rps, status="hit")
                return data
        except requests.RequestException as exc:  # pragma: no cover - network
            last_failure_details = {
                "reason": "network_error",
                "status": None,
                "error": str(exc),
            }
            if attempt >= total_attempts:
                logger.error(
                    "request_error",
                    url=url,
                    error=str(exc),
                    attempt=attempt,
                    total_attempts=total_attempts,
                    rps=cfg.rps,
                )
                logger.debug(
                    "request_fail",
                    url=url,
                    status=None,
                    total_attempts=total_attempts,
                    rps=cfg.rps,
                    **details_excluding("status"),
                )
                return None
            if cfg.delay > 0:
                sleep(cfg.delay)
            continue

    return None


def _ensure_cache(ttl: float) -> TTLCache[str, _CacheEntry]:
    """Return the shared cache instance, recreating it when TTL changes."""

    global _CACHE
    cache = _CACHE
    if cache is None or cache.ttl != ttl:
        cache = _CACHE = TTLCache(maxsize=1024, ttl=ttl)
    return cache


def _store_cache_miss(
    url: str, cfg: PubChemCfg, outcome: str, details: dict[str, Any] | None = None
) -> None:
    """Persist a cached miss outcome for ``url`` including optional details."""

    details_copy = dict(details) if details else None
    with _CACHE_LOCK:
        cache = _ensure_cache(cfg.cache_ttl)
        cache[url] = _CacheEntry(payload=None, outcome=outcome, details=details_copy)
    log_data: dict[str, Any] = {
        "url": url,
        "rps": cfg.rps,
        "status": "miss",
        "outcome": outcome,
    }
    if details_copy:
        log_data.update(details_copy)
    logger.info("cache_set", **log_data)


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
