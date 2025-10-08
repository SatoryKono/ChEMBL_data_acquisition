"""Low-level PubChem client utilities."""

from __future__ import annotations

import json

from dataclasses import dataclass
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from numbers import Real
from threading import Lock
from time import monotonic
from typing import Any, Mapping, cast
from urllib.parse import quote

import requests
from cachetools import TTLCache
from requests import Session

from ..config.models import ApiCfg, PubChemCfg, RetryCfg
from ..config.runtime import session_with_retry
from ..common.log import logger
from ..common.rate_limiter import get_limiter, sleep

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

    base = cfg.base.rstrip("/")
    contains_path_separators = "/" in smiles or "\\" in smiles

    if not contains_path_separators:
        safe_smiles = url_encode(smiles)
        url = f"{base}/compound/smiles/{safe_smiles}/cids/JSON"
        response = make_request(url, cfg)
        if not response:
            return None
    else:
        url = f"{base}/compound/smiles/cids/JSON"
        response = make_request(url, cfg, method="POST", payload={"smiles": smiles})
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


def _build_cache_key(method: str, url: str, data: Mapping[str, Any] | None = None) -> str:
    """Return a stable cache key for ``method``/``url``/``data`` combinations."""

    method_upper = method.upper()
    if not data:
        return f"{method_upper} {url}"
    try:
        serialised = json.dumps(data, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    except TypeError:
        serialised = repr(sorted(data.items()))
    return f"{method_upper} {url}?{serialised}"


def make_request(
    url: str,
    cfg: PubChemCfg,
    *,
    method: str = "GET",
    payload: Mapping[str, Any] | None = None,
) -> dict[str, Any] | None:
    """Make an HTTP request and return parsed JSON."""

    global _CACHE
    method_upper = method.upper()
    cache_key = _build_cache_key(method_upper, url, payload)
    timeout_retry_in: float | None = None
    with _CACHE_LOCK:
        cache = _ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
        cached = cache.get(cache_key) if cache is not None else None
        if (
            cached is not None
            and not cached.is_hit
            and cached.outcome == "timeout"
            and cached.details
        ):
            stored_at = cached.details.get("timeout_stored_at")
            retry_after = cached.details.get("timeout_retry_after")
            if isinstance(stored_at, Real) and isinstance(retry_after, Real):
                elapsed = monotonic() - float(stored_at)
                if elapsed < float(retry_after):
                    timeout_retry_in = float(retry_after) - elapsed
                else:
                    cache.pop(cache_key, None)
                    cached = None
    if cached is not None:
        if cached.is_hit:
            logger.debug("cache_hit", url=url, rps=cfg.rps, status="hit", method=method_upper)
            return cast(dict[str, Any], cached.payload)
        miss_details: dict[str, Any] = {}
        for key, value in (cached.details or {}).items():
            if key == "timeout_stored_at":
                continue
            if key == "status":
                miss_details["http_status"] = value
            else:
                miss_details[key] = value
        if timeout_retry_in is not None:
            miss_details["timeout_retry_in"] = timeout_retry_in
        logger.debug(
            "cache_hit",
            url=url,
            rps=cfg.rps,
            status="miss",
            method=method_upper,
            outcome=cached.outcome,
            **miss_details,
        )
        return None
    logger.debug("cache_miss", url=url, rps=cfg.rps, status="miss", method=method_upper)

    api_cfg = ApiCfg(user_agent=cfg.user_agent)

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
                method=method_upper,
                attempt=attempt,
                total_attempts=total_attempts,
                rps=cfg.rps,
                **(last_failure_details or {}),
            )
            _store_cache_miss(
                cache_key,
                cfg,
                "timeout",
                last_failure_details,
                url=url,
            )
            logger.debug(
                "request_fail",
                url=url,
                status=None,
                total_attempts=total_attempts,
                method=method_upper,
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
            session = get_session(api_cfg)
            request_kwargs: dict[str, Any] = {
                "timeout": (cfg.timeout_connect, cfg.timeout_read),
            }
            if method_upper != "GET" and payload is not None:
                request_kwargs["data"] = payload
            with session.request(method_upper, url, **request_kwargs) as response:
                status = response.status_code
                if status == 404:
                    last_failure_details = {"reason": "not_found", "status": status}
                    logger.info(
                        "request_not_found",
                        url=url,
                        method=method_upper,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    _store_cache_miss(
                        cache_key,
                        cfg,
                        "not_found",
                        last_failure_details,
                        url=url,
                    )
                    logger.debug(
                        "request_fail",
                        url=url,
                        total_attempts=total_attempts,
                        method=method_upper,
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
                    warning_context: dict[str, Any] = {
                        "url": url,
                        "method": method_upper,
                        "status": status,
                        "attempt": attempt,
                        "total_attempts": total_attempts,
                        "rps": cfg.rps,
                    }
                    if retry_after is not None:
                        warning_context["retry_after"] = retry_after
                    logger.warning(
                        event_name,
                        **warning_context,
                    )
                    if attempt >= total_attempts:
                        _store_cache_miss(
                            cache_key,
                            cfg,
                            reason,
                            last_failure_details,
                            url=url,
                        )
                        logger.debug(
                            "request_fail",
                            url=url,
                            total_attempts=total_attempts,
                            method=method_upper,
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
                    if status == 400:
                        last_failure_details = {
                            "reason": "invalid_identifier",
                            "status": status,
                        }
                        logger.info(
                            "request_invalid_identifier",
                            url=url,
                            method=method_upper,
                            rps=cfg.rps,
                            status=status,
                            **details_excluding("status"),
                        )
                        _store_cache_miss(
                            cache_key,
                            cfg,
                            "invalid_identifier",
                            last_failure_details,
                            url=url,
                        )
                    else:
                        last_failure_details = {
                            "reason": "unexpected_status",
                            "status": status,
                        }
                        logger.warning(
                            "request_unexpected_status",
                            url=url,
                            method=method_upper,
                            rps=cfg.rps,
                            status=status,
                            **details_excluding("status"),
                        )
                    logger.debug(
                        "request_fail",
                        url=url,
                        total_attempts=total_attempts,
                        method=method_upper,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    return None

                try:
                    response.raise_for_status()
                    response_data = cast(dict[str, Any], response.json())
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
                            method=method_upper,
                            rps=cfg.rps,
                        )
                        logger.debug(
                            "request_fail",
                            url=url,
                            total_attempts=total_attempts,
                            method=method_upper,
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
                        method=method_upper,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    logger.debug(
                        "request_fail",
                        url=url,
                        total_attempts=total_attempts,
                        method=method_upper,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    return None

                logger.debug(
                    "request_ok",
                    url=url,
                    status=status,
                    method=method_upper,
                    rps=cfg.rps,
                )
                with _CACHE_LOCK:
                    cache = _ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
                    cache[cache_key] = _CacheEntry(payload=response_data, outcome="hit")
                logger.debug(
                    "cache_set",
                    url=url,
                    rps=cfg.rps,
                    status="hit",
                    method=method_upper,
                )
                return response_data
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
                    method=method_upper,
                    rps=cfg.rps,
                )
                logger.debug(
                    "request_fail",
                    url=url,
                    status=None,
                    total_attempts=total_attempts,
                    method=method_upper,
                    rps=cfg.rps,
                    **details_excluding("status"),
                )
                return None
            if cfg.delay > 0:
                sleep(cfg.delay)
            continue

    return None


def _ensure_cache(ttl: float, maxsize: int) -> TTLCache[str, _CacheEntry]:
    """Return the shared cache instance, recreating it when TTL or maxsize change."""

    global _CACHE
    cache = _CACHE
    if cache is None or cache.ttl != ttl or cache.maxsize != maxsize:
        cache = _CACHE = TTLCache(maxsize=maxsize, ttl=ttl)
    return cache


def _store_cache_miss(
    cache_key: str,
    cfg: PubChemCfg,
    outcome: str,
    details: dict[str, Any] | None = None,
    *,
    url: str | None = None,
) -> None:
    """Persist a cached miss outcome for ``url`` including optional details."""

    details_data = dict(details) if details else {}
    cached = False
    log_details: dict[str, Any] = {}
    with _CACHE_LOCK:
        cache = _ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
        if outcome == "timeout":
            base_backoff = cfg.backoff_initial_seconds if cfg.backoff_initial_seconds > 0 else cfg.delay
            retry_after_hint = details_data.get("retry_after")
            if isinstance(retry_after_hint, Real) and float(retry_after_hint) > 0:
                hint_value = float(retry_after_hint)
                if base_backoff and base_backoff > 0:
                    base_backoff = max(base_backoff, hint_value)
                else:
                    base_backoff = hint_value
            existing = cache.get(cache_key)
            if (
                existing is not None
                and not existing.is_hit
                and existing.outcome == "timeout"
                and existing.details
            ):
                previous_retry = existing.details.get("timeout_retry_after")
                if isinstance(previous_retry, Real):
                    doubled = float(previous_retry) * 2
                    if base_backoff and base_backoff > 0:
                        base_backoff = max(base_backoff, doubled)
                    else:
                        base_backoff = doubled
            max_backoff = cfg.timeout_seconds if cfg.timeout_seconds and cfg.timeout_seconds > 0 else None
            effective_backoff = base_backoff if base_backoff and base_backoff > 0 else None
            if effective_backoff is not None:
                if max_backoff is not None:
                    effective_backoff = min(effective_backoff, max_backoff)
                stored_at = monotonic()
                details_data.update(
                    {
                        "timeout": True,
                        "timeout_retry_after": effective_backoff,
                        "timeout_stored_at": stored_at,
                    }
                )
                cache[cache_key] = _CacheEntry(
                    payload=None,
                    outcome=outcome,
                    details=details_data.copy(),
                )
                cached = True
                log_details = {
                    key: value
                    for key, value in details_data.items()
                    if key != "timeout_stored_at"
                }
            else:
                cache.pop(cache_key, None)
                log_details = details_data.copy()
        else:
            cache[cache_key] = _CacheEntry(
                payload=None,
                outcome=outcome,
                details=details_data.copy() if details_data else None,
            )
            cached = True
            log_details = details_data.copy()
    log_data: dict[str, Any] = {
        "url": url or cache_key,
        "rps": cfg.rps,
        "status": "miss",
        "outcome": outcome,
    }
    if log_details:
        log_data.update(log_details)
    event = "cache_set" if cached else "cache_skip"
    logger.debug(event, **log_data)


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
        f"{base}/compound/cid/{validated}/property/MolecularFormula,IUPACName,SMILES,"
        "ConnectivitySMILES,InChI,InChIKey/JSON"
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
        cast(str | None, item.get("SMILES")),
        cast(str | None, item.get("ConnectivitySMILES")),
        cast(str | None, item.get("InChI")),
        cast(str | None, item.get("InChIKey")),
    )
