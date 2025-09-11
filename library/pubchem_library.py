"""PubChem API client utilities.

This module provides functions to interact with the PubChem REST API.
The implementation is a Python translation of a PowerQuery script.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional
from urllib.parse import quote

import requests
from cachetools import LRUCache  # type: ignore[import-untyped]
from requests import Session

from .config import ApiCfg, PubChemCfg, RetryCfg, session_with_retry
from .rate_limiter import get_limiter
from .log import logger

_CACHE: LRUCache[str, Dict[str, Any]] = LRUCache(maxsize=1024)

_session: Session = session_with_retry(ApiCfg(), RetryCfg())


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


def _cids_from_identifier_list(data: Dict[str, Any]) -> List[str]:
    """Extract CIDs from a JSON ``IdentifierList`` structure."""
    return [str(cid) for cid in data.get("IdentifierList", {}).get("CID", [])]


def get_cid_from_smiles(smiles: str, cfg: PubChemCfg) -> Optional[str]:
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


def get_cid_from_inchi(inchi: str, cfg: PubChemCfg) -> Optional[str]:
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


def get_cid_from_inchikey(inchikey: str, cfg: PubChemCfg) -> Optional[str]:
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


def make_request(url: str, cfg: PubChemCfg) -> Optional[Dict[str, Any]]:
    """Make an HTTP GET request and return parsed JSON."""
    if url in _CACHE:
        logger.info("cache_hit", extra={"stage": "cache_hit", "url": url})
        return _CACHE[url]
    logger.info("cache_miss", extra={"stage": "cache_miss", "url": url})

    for attempt in range(1, cfg.retries + 1):
        event = "request_start" if attempt == 1 else "request_retry"
        logger.info(event, extra={"stage": event, "url": url, "attempt": attempt})
        get_limiter("pubchem", cfg.rps, cfg.burst).acquire()
        try:
            response = _session.get(
                url, timeout=(cfg.timeout_connect, cfg.timeout_read)
            )
        except requests.RequestException as exc:  # pragma: no cover - network
            if attempt >= cfg.retries:
                logger.error("HTTP request failed for url %s: %s", url, exc)
                logger.info("request_fail", extra={"stage": "request_fail", "url": url})
                return None
            continue

        if response.status_code in (404, 400):
            logger.warning("Request returned %d for url %s", response.status_code, url)
            logger.info(
                "request_fail",
                extra={
                    "stage": "request_fail",
                    "url": url,
                    "status": response.status_code,
                },
            )
            return None
        try:
            response.raise_for_status()
            data = response.json()
        except requests.RequestException as exc:  # pragma: no cover - network
            if attempt >= cfg.retries:
                logger.error("HTTP request failed for url %s: %s", url, exc)
                logger.info("request_fail", extra={"stage": "request_fail", "url": url})
                return None
            continue
        except ValueError:
            logger.warning("Non-JSON response for url %s", url)
            logger.info(
                "request_fail",
                extra={
                    "stage": "request_fail",
                    "url": url,
                    "status": response.status_code,
                },
            )
            return None

        logger.info(
            "request_ok",
            extra={"stage": "request_ok", "url": url, "status": response.status_code},
        )
        _CACHE[url] = data
        logger.info("cache_set", extra={"stage": "cache_set", "url": url})
        return data
    return None


def validate_cid(cid: str) -> Optional[str]:
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


def _extract_cids(bindings: List[Dict[str, Any]]) -> List[str]:
    """Extract CIDs from API bindings."""
    cids: List[str] = []
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


def get_cid(compound_name: str, cfg: PubChemCfg) -> Optional[str]:
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


def get_all_cid(compound_name: str, cfg: PubChemCfg) -> Optional[str]:
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


def get_standard_name(cid: str, cfg: PubChemCfg) -> Optional[str]:
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
    return info[0].get("Title")


@dataclass
class Properties:
    """Chemical properties for a PubChem compound."""

    IUPACName: str
    MolecularFormula: str
    iSMILES: str
    cSMILES: str
    InChI: str
    InChIKey: str


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
        Chemical property record. Missing values are returned as ``"Not Found"``.
    """
    validated = validate_cid(cid)
    if not validated:
        return Properties(
            "Not Found", "Not Found", "Not Found", "Not Found", "Not Found", "Not Found"
        )
    base = cfg.base.rstrip("/")
    url = (
        f"{base}/compound/cid/{validated}/property/MolecularFormula,IUPACName,IsomericSMILES,"
        "CanonicalSMILES,InChI,InChIKey/JSON"
    )
    response = make_request(url, cfg)
    if not response:
        return Properties(
            "Not Found", "Not Found", "Not Found", "Not Found", "Not Found", "Not Found"
        )
    props = response.get("PropertyTable", {}).get("Properties", [])
    if not props:
        return Properties(
            "Not Found", "Not Found", "Not Found", "Not Found", "Not Found", "Not Found"
        )
    item = props[0]
    return Properties(
        item.get("IUPACName", "Not Found"),
        item.get("MolecularFormula", "Not Found"),
        item.get("IsomericSMILES", "Not Found"),
        item.get("CanonicalSMILES", "Not Found"),
        item.get("InChI", "Not Found"),
        item.get("InChIKey", "Not Found"),
    )


def process_compound(compound_name: str, cfg: PubChemCfg) -> Dict[str, str]:
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
    props = (
        get_properties(cid, cfg)
        if cid
        else Properties(
            "Not Found", "Not Found", "Not Found", "Not Found", "Not Found", "Not Found"
        )
    )
    return {
        "Name": compound_name,
        "CID": cid or "Not Found",
        "Standard Name": standard or "Not Found",
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
    "Properties",
]
