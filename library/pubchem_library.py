"""PubChem API client utilities.

This module provides functions to interact with the PubChem REST API.
The implementation is a Python translation of a PowerQuery script.
"""

from __future__ import annotations

import time
from collections.abc import Mapping, MutableMapping
from dataclasses import dataclass
from typing import Any, cast
from urllib.parse import quote

import requests
from cachetools import TTLCache
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .config import ApiCfg, PubChemCfg, RetryCfg, session_with_retry
from .log import logger
from .rate_limiter import get_limiter, sleep

# Cache is initialised lazily to allow configuration of the TTL via
# :class:`PubChemCfg`. The cache is recreated when the TTL changes.
_CACHE: TTLCache[str, dict[str, Any]] | None = None

# Shared session with placeholder user agent; production code should call
# :func:`init_session` to supply real contact details.
def _build_session(api: ApiCfg, retry: RetryCfg) -> Session:
    """Return HTTP session for PubChem requests disabling automatic retries."""

    session = session_with_retry(api, retry)
    adapter = HTTPAdapter(max_retries=Retry(total=0, allowed_methods=None))
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


_session: Session = _build_session(
    ApiCfg(user_agent="chembl-da/0.1 (mailto:contact@example.org)"), RetryCfg()
)


def init_session(api: ApiCfg, retry: RetryCfg) -> None:
    """Initialise the shared HTTP session.

    Parameters
    ----------
    api:
        Global API settings providing the ``User-Agent`` header.
    retry:
        Retry configuration applied to all requests.
    """

    global _session
    _session = _build_session(api, retry)


def url_encode(text: str) -> str:
    """URL-encode *text* for safe usage in HTTP requests.

    Parameters
    ----------
    text: str
        The string to encode.

    Returns
    -------
    str
        URL-encoded string.

    """
    return quote(text, safe="")


def _cids_from_identifier_list(data: dict[str, Any]) -> list[str]:
    """Extract CIDs from a JSON ``IdentifierList`` structure."""
    return [str(cid) for cid in data.get("IdentifierList", {}).get("CID", [])]


def _normalise_identifier(value: Any, *, uppercase: bool = False) -> str | None:
    """Return ``value`` normalised as a stripped string."""

    if value is None:
        return None
    try:
        if value != value:  # type: ignore[comparison-overlap]
            return None
    except Exception:  # pragma: no cover - defensive
        pass
    text = str(value).strip()
    if not text:
        return None
    return text.upper() if uppercase else text


def _select_primary_cid(
    candidates: list[str],
    *,
    chembl_id: str | None,
    identifier: str,
    value: str | None,
) -> str | None:
    """Return the first unique CID logging alternatives when present."""

    unique: list[str] = []
    seen: set[str] = set()
    for cid in candidates:
        cid_value = cid.strip()
        if not cid_value or cid_value in seen:
            continue
        unique.append(cid_value)
        seen.add(cid_value)
    if not unique:
        return None
    primary = unique[0]
    if len(unique) > 1:
        logger.info(
            "pubchem_multiple_cid",
            chembl_id=chembl_id,
            identifier=identifier,
            value=value,
            cid=primary,
            alternatives=unique[1:],
        )
    return primary


def get_cid_from_smiles(smiles: str, cfg: PubChemCfg) -> str | None:
    """Retrieve PubChem CID(s) for a SMILES string.

    Parameters
    ----------
    smiles: str
        SMILES representation of a compound.
    cfg:
        API configuration providing base URL and timeouts.

    Returns
    -------
    str or None
        Pipe-separated list of CIDs or ``None`` if the structure is
        unknown to PubChem.

    """
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
    """Retrieve PubChem CID(s) for an InChI string.

    Parameters
    ----------
    inchi:
        InChI representation of a compound.
    cfg:
        API configuration providing base URL and timeouts.

    Returns
    -------
    str or None
        Pipe-separated list of CIDs or ``None`` if the structure is
        unknown to PubChem.

    """
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
    """Retrieve PubChem CID(s) for an InChIKey.

    Parameters
    ----------
    inchikey:
        InChIKey representation of a compound.
    cfg:
        API configuration providing base URL and timeouts.

    Returns
    -------
    str or None
        Pipe-separated list of CIDs or ``None`` if the structure is
        unknown to PubChem.
    """
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
    """Make an HTTP GET request and return parsed JSON.

    Parameters
    ----------
    url: str
        Endpoint to query.
    cfg: PubChemCfg
        API configuration. ``cfg.delay`` specifies the pause between retry
        attempts.

    Returns
    -------
    dict[str, Any] or None
        Parsed JSON response on success, otherwise ``None`` when all retries
        are exhausted or a non-recoverable error occurs.
    """
    global _CACHE
    if _CACHE is None or _CACHE.ttl != cfg.cache_ttl:
        # Initialise or refresh the cache with the configured TTL.
        _CACHE = TTLCache(maxsize=1024, ttl=cfg.cache_ttl)

    cached = _CACHE.get(url)
    if cached is not None:
        logger.info("cache_hit", url=url, rps=cfg.rps, status="hit")
        return cast(dict[str, Any], cached)
    logger.info("cache_miss", url=url, rps=cfg.rps, status="miss")

    for attempt in range(1, cfg.retries + 1):
        event = "request_start" if attempt == 1 else "request_retry"
        logger.info(event, url=url, attempt=attempt, rps=cfg.rps)
        get_limiter("pubchem", cfg.rps, cfg.burst).acquire()
        try:
            response = _session.get(
                url, timeout=(cfg.timeout_connect, cfg.timeout_read)
            )
        except requests.RequestException as exc:  # pragma: no cover - network
            if attempt >= cfg.retries:
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
            sleep(cfg.delay)
            continue

        if response.status_code in (404, 400):
            logger.warning(
                "request_unexpected_status",
                url=url,
                status=response.status_code,
                rps=cfg.rps,
            )
            logger.info(
                "request_fail",
                url=url,
                status=response.status_code,
                rps=cfg.rps,
            )
            return None
        try:
            response.raise_for_status()
            data = cast(dict[str, Any], response.json())
        except requests.RequestException as exc:  # pragma: no cover - network
            if attempt >= cfg.retries:
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
                    status=response.status_code,
                    rps=cfg.rps,
                )
                return None
            sleep(cfg.delay)
            continue
        except ValueError:
            logger.warning(
                "response_not_json",
                url=url,
                status=response.status_code,
                rps=cfg.rps,
            )
            logger.info(
                "request_fail",
                url=url,
                status=response.status_code,
                rps=cfg.rps,
            )
            return None

        logger.info(
            "request_ok",
            url=url,
            status=response.status_code,
            rps=cfg.rps,
        )
        assert _CACHE is not None  # for type checker; cache initialised above
        _CACHE[url] = data
        logger.info("cache_set", url=url, rps=cfg.rps)
        return data
    return None


def validate_cid(cid: str) -> str | None:
    """Validate PubChem CID.

    Parameters
    ----------
    cid: str
        Candidate CID.

    Returns
    -------
    str or None
        ``cid`` if valid, otherwise ``None`` when the identifier is empty or
        represents an invalid placeholder (``"0"`` or ``"-1"``).

    """
    if cid in {"", "0", "-1"}:
        return None
    return cid


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
    """Retrieve PubChem CID(s) for *compound_name* (exact match).

    Parameters
    ----------
    compound_name: str
        Compound name to query.
    cfg:
        API configuration providing base URL and timeouts.

    Returns
    -------
    str or None
        Pipe-separated list of CIDs or ``None`` if not found.

    """
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
    """Retrieve PubChem CID(s) for *compound_name* (partial match).

    Parameters
    ----------
    compound_name: str
        Compound name to query.
    cfg:
        API configuration providing base URL and timeouts.

    Returns
    -------
    str or None
        Pipe-separated list of CIDs or ``None`` if not found.
    """
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


def get_standard_name(cid: str, cfg: PubChemCfg) -> str | None:
    """Retrieve the standard compound name for a given CID.

    Parameters
    ----------
    cid: str
        Candidate CID.
    cfg:
        API configuration providing base URL and timeouts.

    Returns
    -------
    str or None
        Standard compound name or ``None`` if ``cid`` is invalid or unknown.
    """
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
    """Retrieve chemical properties for a compound by CID.

    Parameters
    ----------
    cid: str
        Candidate CID.
    cfg:
        API configuration providing base URL and timeouts.

    Returns
    -------
    Properties
        Chemical property record. Missing values are returned as ``None``.
    """
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


@dataclass(frozen=True)
class PubChemRecord:
    """Resolution result combining CID and associated properties."""

    cid: str | None
    properties: Properties
    resolved_by: str | None


def _blank_properties() -> Properties:
    """Return an empty :class:`Properties` instance."""

    return Properties(None, None, None, None, None, None)


def _cache_key(kind: str, value: str) -> str:
    """Return cache key for *kind* and *value*."""

    return f"{kind}:{value}"


def _remaining_time(deadline: float | None) -> float | None:
    """Return seconds remaining until ``deadline`` or ``None`` when unlimited."""

    if deadline is None:
        return None
    return max(0.0, deadline - time.monotonic())


def _backoff_delay(delay: float, cfg: PubChemCfg) -> float:
    """Return the next backoff interval."""

    if delay > 0:
        return delay
    if cfg.backoff_initial_seconds > 0:
        return cfg.backoff_initial_seconds
    if cfg.delay > 0:
        return cfg.delay
    return 0.5


def _record_from_cache(
    cache: MutableMapping[str, PubChemRecord] | None, kind: str, value: str
) -> PubChemRecord | None:
    if cache is None:
        return None
    return cache.get(_cache_key(kind, value))


def _store_record(
    cache: MutableMapping[str, PubChemRecord] | None,
    kind: str,
    value: str,
    record: PubChemRecord,
) -> None:
    if cache is None:
        return
    cache[_cache_key(kind, value)] = record


def _properties_for_cid(
    cache: MutableMapping[str, PubChemRecord] | None, cid: str, cfg: PubChemCfg
) -> Properties:
    cached = _record_from_cache(cache, "cid", cid)
    if cached and cached.cid:
        return cached.properties
    return get_properties(cid, cfg)


def _identifier_urls(kind: str, value: str, base: str) -> list[str]:
    encoded = url_encode(value)
    if kind == "smiles":
        return [f"{base}/compound/smiles/{encoded}/cids/JSON"]
    if kind == "inchikey":
        return [f"{base}/compound/inchikey/{encoded}/cids/JSON"]
    if kind == "inchi":
        return [f"{base}/compound/inchi/{encoded}/cids/JSON"]
    if kind == "pref_name":
        return [
            f"{base}/compound/name/{encoded}/cids/JSON",
            f"{base}/compound/name/{encoded}/cids/JSON?name_type=word",
        ]
    raise ValueError(f"unsupported resolution step: {kind}")


def _request_with_backoff(
    url: str, cfg: PubChemCfg, *, deadline: float | None
) -> dict[str, Any] | None:
    """Return JSON response handling backoff for rate limiting and 5xx errors."""

    global _CACHE
    if _CACHE is None or _CACHE.ttl != cfg.cache_ttl:
        _CACHE = TTLCache(maxsize=1024, ttl=cfg.cache_ttl)

    cached = _CACHE.get(url)
    if cached is not None:
        logger.info("cache_hit", url=url, rps=cfg.rps, status="hit")
        return cast(dict[str, Any], cached)
    logger.info("cache_miss", url=url, rps=cfg.rps, status="miss")

    delay = cfg.backoff_initial_seconds
    attempt = 0
    while True:
        attempt += 1
        event = "request_start" if attempt == 1 else "request_retry"
        logger.info(event, url=url, attempt=attempt, rps=cfg.rps)
        get_limiter("pubchem", cfg.rps, cfg.burst).acquire()
        try:
            response = _session.get(
                url, timeout=(cfg.timeout_connect, cfg.timeout_read)
            )
        except requests.RequestException as exc:  # pragma: no cover - network
            remaining = _remaining_time(deadline)
            logger.warning(
                "request_error",
                url=url,
                error=str(exc),
                attempt=attempt,
                rps=cfg.rps,
            )
            if remaining is not None and remaining <= 0:
                logger.info("request_fail", url=url, status=None, rps=cfg.rps)
                return None
            wait = _backoff_delay(delay, cfg)
            if remaining is not None:
                wait = min(wait, remaining)
            if wait > 0:
                logger.info(
                    "pubchem_retry",
                    url=url,
                    status=None,
                    attempt=attempt,
                    delay=wait,
                    rps=cfg.rps,
                )
                sleep(wait)
            delay = wait * 2 if wait > 0 else delay
            continue

        status = response.status_code
        if status == 404:
            logger.info(
                "request_not_found", url=url, status=status, attempt=attempt, rps=cfg.rps
            )
            logger.info("request_fail", url=url, status=status, rps=cfg.rps)
            return None
        if status == 429 or 500 <= status < 600:
            remaining = _remaining_time(deadline)
            if remaining is not None and remaining <= 0:
                logger.warning(
                    "pubchem_timeout",
                    url=url,
                    status=status,
                    attempt=attempt,
                    rps=cfg.rps,
                )
                logger.info("request_fail", url=url, status=status, rps=cfg.rps)
                return None
            wait = _backoff_delay(delay, cfg)
            if remaining is not None:
                wait = min(wait, remaining)
            if wait > 0:
                logger.info(
                    "pubchem_retry",
                    url=url,
                    status=status,
                    attempt=attempt,
                    delay=wait,
                    rps=cfg.rps,
                )
                sleep(wait)
            delay = wait * 2 if wait > 0 else delay
            continue
        try:
            response.raise_for_status()
            data = cast(dict[str, Any], response.json())
        except requests.HTTPError as exc:  # pragma: no cover - network
            logger.error(
                "request_error",
                url=url,
                error=str(exc),
                attempt=attempt,
                status=status,
                rps=cfg.rps,
            )
            logger.info("request_fail", url=url, status=status, rps=cfg.rps)
            return None
        except ValueError:
            logger.warning(
                "response_not_json",
                url=url,
                status=status,
                rps=cfg.rps,
            )
            logger.info("request_fail", url=url, status=status, rps=cfg.rps)
            return None

        logger.info("request_ok", url=url, status=status, rps=cfg.rps)
        assert _CACHE is not None
        _CACHE[url] = data
        logger.info("cache_set", url=url, rps=cfg.rps)
        return data


def _resolve_identifier(
    kind: str,
    value: str,
    cfg: PubChemCfg,
    *,
    deadline: float | None,
    chembl_id: str | None,
) -> str | None:
    base = cfg.base.rstrip("/")
    for url in _identifier_urls(kind, value, base):
        if deadline is not None and _remaining_time(deadline) == 0:
            return None
        data = _request_with_backoff(url, cfg, deadline=deadline)
        if not data:
            continue
        cid = _select_primary_cid(
            _cids_from_identifier_list(data),
            chembl_id=chembl_id,
            identifier=kind,
            value=value,
        )
        if cid:
            return cid
    return None


def _is_truthy_flag(value: Any) -> bool:
    text = _normalise_identifier(value) or ""
    return text.lower() in {"1", "true", "yes", "y", "t"}


def resolve_pubchem_record(
    row: Mapping[str, Any],
    cfg: PubChemCfg,
    *,
    cid_cache: MutableMapping[str, str | None] | None = None,
    resolution_cache: MutableMapping[str, PubChemRecord] | None = None,
) -> PubChemRecord:
    """Resolve PubChem information for ``row`` according to ``cfg``."""

    if not cfg.enable:
        return PubChemRecord(None, _blank_properties(), None)

    chembl_id = _normalise_identifier(row.get("molecule_chembl_id"), uppercase=True)
    parent_id = _normalise_identifier(
        row.get("parent_molecule_chembl_id"), uppercase=True
    )

    if not cfg.allow_polymer and _is_truthy_flag(row.get("polymer_flag")):
        if cid_cache is not None and chembl_id:
            cid_cache.setdefault(chembl_id, None)
        return PubChemRecord(None, _blank_properties(), None)

    resolution_cache = resolution_cache if resolution_cache is not None else {}

    row_cid = _normalise_identifier(row.get("pubchem_cid"))
    if row_cid:
        if cid_cache is not None and chembl_id:
            cid_cache[chembl_id] = row_cid
        record = _record_from_cache(resolution_cache, "cid", row_cid)
        if record and record.cid:
            return PubChemRecord(record.cid, record.properties, "cached_cid")
        props = _properties_for_cid(resolution_cache, row_cid, cfg)
        record = PubChemRecord(row_cid, props, "cached_cid")
        _store_record(resolution_cache, "cid", row_cid, record)
        return record

    if chembl_id and cid_cache is not None and chembl_id in cid_cache:
        cached_cid = cid_cache[chembl_id]
        if cached_cid:
            record = _record_from_cache(resolution_cache, "cid", cached_cid)
            if record and record.cid:
                return PubChemRecord(record.cid, record.properties, "cached_cid")
            props = _properties_for_cid(resolution_cache, cached_cid, cfg)
            record = PubChemRecord(cached_cid, props, "cached_cid")
            _store_record(resolution_cache, "cid", cached_cid, record)
            return record
        return PubChemRecord(None, _blank_properties(), None)

    if cfg.use_parent_for_salts and parent_id and cid_cache is not None:
        parent_cid = cid_cache.get(parent_id)
        if parent_cid:
            record = _record_from_cache(resolution_cache, "cid", parent_cid)
            if record and record.cid:
                if chembl_id:
                    cid_cache[chembl_id] = parent_cid
                return PubChemRecord(record.cid, record.properties, "cached_cid")
            props = _properties_for_cid(resolution_cache, parent_cid, cfg)
            record = PubChemRecord(parent_cid, props, "cached_cid")
            _store_record(resolution_cache, "cid", parent_cid, record)
            if chembl_id:
                cid_cache[chembl_id] = parent_cid
            return record

    deadline = (
        time.monotonic() + cfg.timeout_seconds if cfg.timeout_seconds > 0 else None
    )

    for step in cfg.resolve_order:
        if step == "cached_cid":
            continue
        if step == "smiles":
            candidate = _normalise_identifier(row.get("canonical_smiles"))
        elif step == "inchikey":
            candidate = _normalise_identifier(
                row.get("standard_inchi_key"), uppercase=True
            )
        elif step == "inchi":
            candidate = _normalise_identifier(row.get("standard_inchi"))
        elif step == "pref_name":
            candidate = _normalise_identifier(row.get("pref_name"))
        else:
            candidate = None
        if not candidate:
            continue

        cached_record = _record_from_cache(resolution_cache, step, candidate)
        if cached_record is not None:
            if cached_record.cid and cid_cache is not None and chembl_id:
                cid_cache[chembl_id] = cached_record.cid
            resolved_by = step if cached_record.cid else None
            return PubChemRecord(
                cached_record.cid, cached_record.properties, resolved_by
            )

        cid = _resolve_identifier(
            step, candidate, cfg, deadline=deadline, chembl_id=chembl_id
        )
        if not cid:
            _store_record(
                resolution_cache,
                step,
                candidate,
                PubChemRecord(None, _blank_properties(), None),
            )
            continue

        props = _properties_for_cid(resolution_cache, cid, cfg)
        record = PubChemRecord(cid, props, step)
        _store_record(resolution_cache, step, candidate, record)
        _store_record(resolution_cache, "cid", cid, record)
        if cid_cache is not None and chembl_id:
            cid_cache[chembl_id] = cid
        if cfg.use_parent_for_salts and parent_id and cid_cache is not None:
            cid_cache.setdefault(parent_id, cid)
        return record

    if cid_cache is not None and chembl_id:
        cid_cache.setdefault(chembl_id, None)
    return PubChemRecord(None, _blank_properties(), None)


def process_compound(compound_name: str, cfg: PubChemCfg) -> dict[str, str | None]:
    """Process *compound_name* into a structured record.

    Parameters
    ----------
    compound_name: str
        Name of the compound to look up.
    cfg:
        API configuration providing base URL and timeouts.

    Returns
    -------
    dict
        Dictionary containing compound details.

    """
    cid = get_cid(compound_name, cfg)
    standard = get_standard_name(cid, cfg) if cid else None
    props = get_properties(cid, cfg) if cid else Properties(None, None, None, None, None, None)
    return {
        "Name": compound_name,
        "CID": cid,
        "Standard Name": standard,
        "IUPACName": props.IUPACName,
        "MolecularFormula": props.MolecularFormula,
        "iSMILES": props.iSMILES,
        "cSMILES": props.cSMILES,
        "InChI": props.InChI,
        "InChIKey": props.InChIKey,
    }


__all__ = [
    "url_encode",
    "init_session",
    "make_request",
    "validate_cid",
    "get_cid",
    "get_all_cid",
    "get_standard_name",
    "get_properties",
    "resolve_pubchem_record",
    "process_compound",
    "Properties",
    "PubChemRecord",
]
