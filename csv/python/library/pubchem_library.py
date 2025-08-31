"""PubChem API client utilities.

This module provides functions to interact with the PubChem REST API.
The implementation is a Python translation of a PowerQuery script.
"""
from __future__ import annotations

from dataclasses import dataclass
import logging
import time
from typing import Any, Dict, List, Optional
from urllib.parse import quote

import requests
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


logger = logging.getLogger(__name__)

# A single shared session with retry/backoff for all HTTP calls.  PubChem
# enforces fairly strict rate limits; the retry configuration helps to recover
# from transient failures such as HTTP 5xx responses.
_retry = Retry(
    total=3,
    backoff_factor=1.0,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET"],
)
_session: Session = requests.Session()
_session.mount("http://", HTTPAdapter(max_retries=_retry))
_session.mount("https://", HTTPAdapter(max_retries=_retry))


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


def get_cid_from_smiles(smiles: str) -> Optional[str]:
    """Retrieve PubChem CID(s) for a SMILES string.

    Parameters
    ----------
    smiles: str
        SMILES representation of a compound.

    Returns
    -------
    str or None
        Pipe-separated list of CIDs or ``None`` if the structure is
        unknown to PubChem.
    """

    safe_smiles = url_encode(smiles)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
        f"{safe_smiles}/cids/JSON"
    )
    response = make_request(url)
    if not response:
        return None
    cids = _cids_from_identifier_list(response)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_cid_from_inchi(inchi: str) -> Optional[str]:
    """Retrieve PubChem CID(s) for an InChI string.

    Parameters
    ----------
    inchi:
        InChI representation of a compound.

    Returns
    -------
    str or None
        Pipe-separated list of CIDs or ``None`` if the structure is
        unknown to PubChem.
    """

    safe_inchi = url_encode(inchi)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/"
        f"{safe_inchi}/cids/JSON"
    )
    response = make_request(url)
    if not response:
        return None
    cids = _cids_from_identifier_list(response)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_cid_from_inchikey(inchikey: str) -> Optional[str]:
    """Retrieve PubChem CID(s) for an InChIKey."""

    safe_inchikey = url_encode(inchikey)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
        f"{safe_inchikey}/cids/JSON"
    )
    response = make_request(url)
    if not response:
        return None
    cids = _cids_from_identifier_list(response)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def make_request(url: str, delay: float = 3.0) -> Optional[Dict[str, Any]]:
    """Make an HTTP GET request and return parsed JSON.

    The function sleeps for ``delay`` seconds before issuing the request in
    order to respect PubChem rate limits.  A shared session configured with
    retries is used to automatically retry transient failures.

    Parameters
    ----------
    url:
        Endpoint URL to query.
    delay:
        Time in seconds to wait before making the request. Defaults to three
        seconds.

    Returns
    -------
    dict or None
        Parsed JSON response or ``None`` when the request fails, the server
        returns a non-success status code, or the payload cannot be decoded.
    """

    time.sleep(delay)
    try:
        response = _session.get(url, timeout=10)
        if response.status_code == 404:
            logger.warning("Request returned 404 for url %s", url)
            return None
        if response.status_code == 400:
            logger.warning("Request returned 400 for url %s", url)
            return None
        response.raise_for_status()
        try:
            return response.json()
        except ValueError:
            logger.warning("Non-JSON response for url %s", url)
            return None
    except requests.RequestException as exc:  # pragma: no cover - network
        logger.error("HTTP request failed for url %s: %s", url, exc)
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


def get_cid(compound_name: str) -> Optional[str]:
    """Retrieve PubChem CID(s) for *compound_name* (exact match).

    Parameters
    ----------
    compound_name: str
        Compound name to query.

    Returns
    -------
    str or None
        Pipe-separated list of CIDs or ``None`` if not found.
    """
    safe_name = url_encode(compound_name)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/rdf/query?graph=synonym&"
        f"return=cid&format=json&name={safe_name}"
    )
    response = make_request(url)
    if not response:
        return None
    bindings = response.get("results", {}).get("bindings", [])
    cids = _extract_cids(bindings)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_all_cid(compound_name: str) -> Optional[str]:
    """Retrieve PubChem CID(s) for *compound_name* (partial match)."""
    safe_name = url_encode(compound_name)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/rdf/query?graph=synonym&"
        f"return=cid&format=json&name={safe_name}&contain=true"
    )
    response = make_request(url)
    if not response:
        return None
    bindings = response.get("results", {}).get("bindings", [])
    cids = _extract_cids(bindings)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_standard_name(cid: str) -> Optional[str]:
    """Retrieve the standard compound name for a given CID."""
    validated = validate_cid(cid)
    if not validated:
        return None
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{validated}/description/JSON"
    )
    response = make_request(url)
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


def get_properties(cid: str) -> Properties:
    """Retrieve chemical properties for a compound by CID."""
    validated = validate_cid(cid)
    if not validated:
        return Properties(
            "Not Found", "Not Found", "Not Found", "Not Found", "Not Found", "Not Found"
        )
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{validated}/property/MolecularFormula,IUPACName,IsomericSMILES,"
        "CanonicalSMILES,InChI,InChIKey/JSON"
    )
    response = make_request(url)
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


def process_compound(compound_name: str) -> Dict[str, str]:
    """Process *compound_name* into a structured record.

    Parameters
    ----------
    compound_name: str
        Name of the compound to look up.

    Returns
    -------
    dict
        Dictionary containing compound details.
    """
    cid = get_cid(compound_name)
    standard = get_standard_name(cid) if cid else None
    props = get_properties(cid) if cid else Properties(
        "Not Found", "Not Found", "Not Found", "Not Found", "Not Found", "Not Found"
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
    "make_request",
    "validate_cid",
    "get_cid",
    "get_all_cid",
    "get_standard_name",
    "get_properties",
    "process_compound",
    "Properties",
]
