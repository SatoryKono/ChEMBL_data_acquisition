"""PubChem API client utilities.

This module provides functions to interact with the PubChem REST API.
The implementation is a Python translation of a PowerQuery script.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, cast, Literal
from urllib.parse import quote

import requests
from cachetools import TTLCache
from requests import Session

from .config import ApiCfg, PubChemCfg, RetryCfg, session_with_retry
from .log import logger
from .rate_limiter import get_limiter, sleep

# Cache is initialised lazily to allow configuration of the TTL via
# :class:`PubChemCfg`. The cache is recreated when the TTL changes.
_CACHE: TTLCache[str, dict[str, Any]] | None = None

# Shared session with placeholder user agent; production code should call
# :func:`init_session` to supply real contact details.
_session: Session = session_with_retry(
    ApiCfg(user_agent="chembl-da/0.1 (mailto:contact@example.org)"), RetryCfg()
)

_RESOLVE_KEYS: frozenset[str] = frozenset({"cache", "smiles", "inchikey", "inchi", "pref_name"})


@dataclass(frozen=True)
class PubChemResolution:
    """Result of :func:`resolve_pubchem_record`."""

    cid: str | None
    source: str | None
    status: Literal["resolved", "not_found", "error", "skipped"]
    status_code: int | None
    attempts: tuple[str, ...]


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


def _primary_cid(value: str | None) -> str | None:
    """Return the first CID from a pipe-separated string."""

    if not value:
        return None
    cids = [part.strip() for part in str(value).split("|") if part.strip()]
    return cids[0] if cids else None


def _request_json(url: str, cfg: PubChemCfg, session: Session) -> tuple[dict[str, Any] | None, int | None]:
    """Perform a JSON request handling retryable PubChem status codes."""

    attempts = max(1, cfg.retries)
    backoff = cfg.backoff_initial_seconds
    status_code: int | None = None
    for attempt in range(1, attempts + 1):
        get_limiter("pubchem", cfg.rps, cfg.burst).acquire()
        try:
            response = session.get(url, timeout=cfg.timeout_seconds)
        except requests.RequestException as exc:
            logger.warning(
                "pubchem_request_error",
                url=url,
                attempt=attempt,
                error=str(exc),
            )
            status_code = None
        else:
            status_code = response.status_code
            if status_code == 200:
                try:
                    return cast(dict[str, Any], response.json()), status_code
                except ValueError:
                    logger.warning(
                        "pubchem_response_not_json",
                        url=url,
                        status=status_code,
                    )
                    return None, status_code
            if status_code == 404:
                logger.info("pubchem_not_found", url=url, status=status_code)
                return None, status_code
            if status_code == 429 or 500 <= status_code < 600:
                logger.warning(
                    "pubchem_retryable_status",
                    url=url,
                    status=status_code,
                    attempt=attempt,
                )
            elif 400 <= status_code < 500:
                logger.warning(
                    "pubchem_client_error",
                    url=url,
                    status=status_code,
                )
                return None, status_code
            else:
                logger.warning(
                    "pubchem_unexpected_status",
                    url=url,
                    status=status_code,
                )
                return None, status_code

        if attempt >= attempts:
            break
        if status_code in (429,) or (status_code is not None and status_code >= 500) or status_code is None:
            delay = backoff if backoff > 0 else cfg.delay
            if delay > 0:
                sleep(delay)
            if backoff > 0:
                backoff *= 2
            continue
        break
    return None, status_code


def _resolve_identifier(
    kind: str,
    value: str,
    cfg: PubChemCfg,
    session: Session,
) -> tuple[list[str], int | None]:
    """Resolve *value* of type *kind* into PubChem CIDs."""

    base = cfg.base.rstrip("/")
    if kind == "smiles":
        url = f"{base}/compound/smiles/{url_encode(value)}/cids/JSON"
    elif kind == "inchikey":
        url = f"{base}/compound/inchikey/{url_encode(value)}/cids/JSON"
    elif kind == "inchi":
        url = f"{base}/compound/inchi/{url_encode(value)}/cids/JSON"
    else:  # pragma: no cover - defensive programming
        raise ValueError(f"Unsupported identifier kind: {kind}")
    data, status = _request_json(url, cfg, session)
    if not data:
        return [], status
    cids = sorted(set(_cids_from_identifier_list(data)))
    return cids, status


def _resolve_pref_name(
    value: str,
    cfg: PubChemCfg,
    session: Session,
) -> tuple[list[str], int | None]:
    """Resolve preferred names using exact and partial synonym searches."""

    rdf_base = cfg.base.rstrip("/").rsplit("/", 1)[0] + "/rdf"
    encoded = url_encode(value)
    urls = [
        f"{rdf_base}/query?graph=synonym&return=cid&format=json&name={encoded}",
        f"{rdf_base}/query?graph=synonym&return=cid&format=json&name={encoded}&contain=true",
    ]
    last_status: int | None = None
    for index, url in enumerate(urls):
        data, status = _request_json(url, cfg, session)
        last_status = status
        if data:
            bindings = data.get("results", {}).get("bindings", [])
            cids = sorted(set(_extract_cids(bindings)))
            if cids:
                return cids, status
        if status not in (404, 200):
            return [], status
        if index == 0 and status in (200, 404):
            continue
        break
    return [], last_status


def resolve_pubchem_record(
    identifiers: Mapping[str, str | None],
    cfg: PubChemCfg,
    *,
    session: Session | None = None,
) -> PubChemResolution:
    """Resolve a PubChem record using configured fallbacks."""

    if not cfg.enable:
        return PubChemResolution(None, None, "skipped", None, tuple())

    molecule_type = (identifiers.get("molecule_type") or "").strip().lower()
    if molecule_type == "polymer" and not cfg.allow_polymer:
        return PubChemResolution(None, None, "skipped", None, tuple())

    session = session or _session
    attempts: list[str] = []
    last_status: int | None = None
    for source in cfg.resolve_order:
        if source not in _RESOLVE_KEYS:
            raise ValueError(f"Unsupported resolve step: {source}")
        attempts.append(source)
        if source == "cache":
            cached = identifiers.get("cache")
            cid = _primary_cid(cached) if cached is not None else None
            if cid:
                return PubChemResolution(cid, source, "resolved", None, tuple(attempts))
            continue
        value = identifiers.get(source)
        if not value:
            continue
        if source == "pref_name":
            cids, status = _resolve_pref_name(value, cfg, session)
        else:
            cids, status = _resolve_identifier(source, value, cfg, session)
        last_status = status
        if cids:
            return PubChemResolution(cids[0], source, "resolved", status, tuple(attempts))
        if status in (404, 200):
            continue
        return PubChemResolution(None, source, "error", status, tuple(attempts))

    return PubChemResolution(None, None, "not_found", last_status, tuple(attempts))


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
    "process_compound",
    "resolve_pubchem_record",
    "PubChemResolution",
    "Properties",
]
