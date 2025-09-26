"""PubChem API client utilities.

This module provides functions to interact with the PubChem REST API.
The implementation is a Python translation of a PowerQuery script.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Callable, Mapping, cast
from urllib.parse import quote

import requests
from cachetools import TTLCache
from requests import Session

from .config import ApiCfg, PubChemCfg, RetryCfg, session_with_retry
from .log import logger
from .rate_limiter import get_limiter, sleep

PUBCHEM_FIELDS = (
    "pubchem_cid",
    "pubchem_iupac_name",
    "pubchem_molecular_formula",
    "pubchem_isomeric_smiles",
    "pubchem_canonical_smiles",
    "pubchem_inchi",
    "pubchem_inchikey",
)

# Cache is initialised lazily to allow configuration of the TTL via
# :class:`PubChemCfg`. The cache is recreated when the TTL changes.
_CACHE: TTLCache[str, dict[str, Any]] | None = None

# Shared session with placeholder user agent; production code should call
# :func:`init_session` to supply real contact details.
_session: Session = session_with_retry(
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
    _session = session_with_retry(api, retry)


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


def _fetch_cids(
    url: str,
    cfg: PubChemCfg,
    *,
    deadline: float | None = None,
) -> tuple[int | None, list[str]]:
    status, payload = make_request(url, cfg, deadline=deadline, return_status=True)
    if status != 200 or payload is None:
        return status, []
    cids = _cids_from_identifier_list(payload)
    return status, sorted(set(cids))


def _fetch_cids_for_name(
    name: str,
    cfg: PubChemCfg,
    *,
    deadline: float | None = None,
) -> tuple[int | None, list[str]]:
    safe_name = url_encode(name)
    rdf_base = cfg.base.rstrip("/").rsplit("/", 1)[0] + "/rdf"
    url = f"{rdf_base}/query?graph=synonym&return=cid&format=json&name={safe_name}"
    status, payload = make_request(url, cfg, deadline=deadline, return_status=True)
    if status != 200 or payload is None:
        return status, []
    bindings = payload.get("results", {}).get("bindings", [])
    cids = _extract_cids(bindings)
    return status, sorted(set(cids))


def get_cid_from_smiles(
    smiles: str, cfg: PubChemCfg, *, deadline: float | None = None
) -> str | None:
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
    status, cids = _fetch_cids(url, cfg, deadline=deadline)
    if status != 200 or not cids:
        return None
    return "|".join(cids)


def get_cid_from_inchi(
    inchi: str, cfg: PubChemCfg, *, deadline: float | None = None
) -> str | None:
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
    status, cids = _fetch_cids(url, cfg, deadline=deadline)
    if status != 200 or not cids:
        return None
    return "|".join(cids)


def get_cid_from_inchikey(
    inchikey: str, cfg: PubChemCfg, *, deadline: float | None = None
) -> str | None:
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
    status, cids = _fetch_cids(url, cfg, deadline=deadline)
    if status != 200 or not cids:
        return None
    return "|".join(cids)


def make_request(
    url: str,
    cfg: PubChemCfg,
    *,
    deadline: float | None = None,
    return_status: bool = False,
) -> dict[str, Any] | tuple[int | None, dict[str, Any] | None] | None:
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
        _CACHE = TTLCache(maxsize=1024, ttl=cfg.cache_ttl)

    cached = _CACHE.get(url)
    if cached is not None:
        logger.info("cache_hit", url=url, rps=cfg.rps, status="hit")
        data = cast(dict[str, Any], cached)
        if return_status:
            return 200, data
        return data
    logger.info("cache_miss", url=url, rps=cfg.rps, status="miss")

    effective_deadline = deadline
    if effective_deadline is None and cfg.timeout_seconds > 0:
        effective_deadline = time.monotonic() + cfg.timeout_seconds

    backoff_delay = cfg.backoff_initial_seconds if cfg.backoff_initial_seconds else cfg.delay
    last_status: int | None = None

    for attempt in range(1, cfg.retries + 1):
        if effective_deadline is not None and time.monotonic() > effective_deadline:
            logger.info("request_timeout", url=url, attempt=attempt, rps=cfg.rps)
            if return_status:
                return 408, None
            return None

        event = "request_start" if attempt == 1 else "request_retry"
        logger.info(event, url=url, attempt=attempt, rps=cfg.rps)
        get_limiter("pubchem", cfg.rps, cfg.burst).acquire()
        try:
            response = _session.get(
                url, timeout=(cfg.timeout_connect, cfg.timeout_read)
            )
        except requests.RequestException as exc:  # pragma: no cover - network
            last_status = None
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
                if return_status:
                    return None, None
                return None
            sleep(cfg.delay)
            continue

        last_status = response.status_code
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
            if return_status:
                return response.status_code, None
            return None

        if response.status_code == 429 or 500 <= response.status_code < 600:
            if attempt >= cfg.retries:
                logger.info(
                    "request_fail",
                    url=url,
                    status=response.status_code,
                    rps=cfg.rps,
                )
                if return_status:
                    return response.status_code, None
                return None
            wait = max(backoff_delay, 0.0)
            if wait:
                sleep(wait)
            if cfg.backoff_initial_seconds:
                backoff_delay = backoff_delay * 2 if backoff_delay else cfg.backoff_initial_seconds
            continue

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
                if return_status:
                    return response.status_code, None
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
            if return_status:
                return response.status_code, None
            return None

        logger.info(
            "request_ok",
            url=url,
            status=response.status_code,
            rps=cfg.rps,
        )
        assert _CACHE is not None
        _CACHE[url] = data
        logger.info("cache_set", url=url, rps=cfg.rps)
        if return_status:
            return response.status_code, data
        return data

    if return_status:
        return last_status, None
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
    if not cid.isdigit():
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


def get_cid(
    compound_name: str, cfg: PubChemCfg, *, deadline: float | None = None
) -> str | None:
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
    response = make_request(url, cfg, deadline=deadline)
    if not response:
        return None
    bindings = response.get("results", {}).get("bindings", [])
    cids = _extract_cids(bindings)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_all_cid(
    compound_name: str, cfg: PubChemCfg, *, deadline: float | None = None
) -> str | None:
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
    response = make_request(url, cfg, deadline=deadline)
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


def _empty_record(write_literal: bool) -> dict[str, str | None]:
    record = {field: None for field in PUBCHEM_FIELDS}
    if write_literal:
        record["pubchem_cid"] = "NOT_FOUND"
    return record


def _build_record(cid: str, props: Properties) -> dict[str, str | None]:
    return {
        "pubchem_cid": cid,
        "pubchem_iupac_name": props.IUPACName,
        "pubchem_molecular_formula": props.MolecularFormula,
        "pubchem_isomeric_smiles": props.iSMILES,
        "pubchem_canonical_smiles": props.cSMILES,
        "pubchem_inchi": props.InChI,
        "pubchem_inchikey": props.InChIKey,
    }


def _normalize_identifier(value: str | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def get_properties(
    cid: str, cfg: PubChemCfg, *, deadline: float | None = None
) -> Properties:
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
    response = make_request(url, cfg, deadline=deadline)
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


def resolve_pubchem_record(
    identifiers: Mapping[str, str | None],
    cfg: PubChemCfg,
    *,
    cache: dict[tuple[str, str], tuple[bool, dict[str, str | None] | None]] | None = None,
) -> dict[str, str | None]:
    """Resolve PubChem information using identifier fallbacks."""

    deadline: float | None = None
    if cfg.timeout_seconds > 0:
        deadline = time.monotonic() + cfg.timeout_seconds

    def cache_get(key: tuple[str, str]) -> tuple[bool, dict[str, str | None] | None] | None:
        if cache is None:
            return None
        return cache.get(key)

    def cache_set(
        key: tuple[str, str],
        found: bool,
        record: dict[str, str | None] | None = None,
    ) -> None:
        if cache is None:
            return
        cache[key] = (found, dict(record) if record is not None else None)

    fetchers: dict[str, Callable[[str], tuple[int | None, list[str]]]] = {
        "smiles": lambda value: _fetch_cids(
            f"{cfg.base.rstrip('/')}/compound/smiles/{url_encode(value)}/cids/JSON",
            cfg,
            deadline=deadline,
        ),
        "inchikey": lambda value: _fetch_cids(
            f"{cfg.base.rstrip('/')}/compound/inchikey/{url_encode(value)}/cids/JSON",
            cfg,
            deadline=deadline,
        ),
        "inchi": lambda value: _fetch_cids(
            f"{cfg.base.rstrip('/')}/compound/inchi/{url_encode(value)}/cids/JSON",
            cfg,
            deadline=deadline,
        ),
        "pref_name": lambda value: _fetch_cids_for_name(
            value,
            cfg,
            deadline=deadline,
        ),
    }

    for method in cfg.resolve_order:
        raw_value = identifiers.get(method)
        value = _normalize_identifier(raw_value)
        if not value:
            continue
        key = (method, value)
        cached = cache_get(key)
        if cached is not None:
            found, record = cached
            if found and record is not None:
                return dict(record)
            if not found:
                continue

        if method == "cid":
            cid = validate_cid(value)
            if not cid:
                cache_set(key, False)
                continue
            props = get_properties(cid, cfg, deadline=deadline)
            record = _build_record(cid, props)
            cache_set(key, True, record)
            cache_set(("cid", cid), True, record)
            return dict(record)

        fetcher = fetchers.get(method)
        if fetcher is None:
            continue
        status, cids = fetcher(value)
        if status == 404 or (status == 200 and not cids):
            cache_set(key, False)
            continue
        if not cids:
            continue
        cid_candidate = validate_cid(cids[0])
        if not cid_candidate:
            cache_set(key, False)
            continue
        props = get_properties(cid_candidate, cfg, deadline=deadline)
        record = _build_record(cid_candidate, props)
        cache_set(key, True, record)
        cache_set(("cid", cid_candidate), True, record)
        return dict(record)

    return _empty_record(cfg.write_not_found_literal)


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
]
