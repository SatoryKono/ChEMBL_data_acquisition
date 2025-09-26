"""PubChem API client utilities.

This module provides functions to interact with the PubChem REST API.
The implementation is a Python translation of a PowerQuery script.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Callable, Mapping, MutableMapping, cast
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

_PROPERTY_FIELDS = (
    "MolecularFormula,IUPACName,IsomericSMILES,CanonicalSMILES,InChI,InChIKey"
)

_IDENTIFIER_PATHS: dict[str, Callable[[str, str], str]] = {
    "cid": lambda base, value: (
        f"{base}/compound/cid/{value}/property/{_PROPERTY_FIELDS}/JSON"
    ),
    "smiles": lambda base, value: (
        f"{base}/compound/smiles/{url_encode(value)}/property/{_PROPERTY_FIELDS}/JSON"
    ),
    "inchikey": lambda base, value: (
        f"{base}/compound/inchikey/{url_encode(value)}/property/{_PROPERTY_FIELDS}/JSON"
    ),
    "inchi": lambda base, value: (
        f"{base}/compound/inchi/{url_encode(value)}/property/{_PROPERTY_FIELDS}/JSON"
    ),
    "pref_name": lambda base, value: (
        f"{base}/compound/name/{url_encode(value)}/property/{_PROPERTY_FIELDS}/JSON"
    ),
}

_NOT_FOUND_LITERAL = "Not Found"


def _missing_literal(cfg: PubChemCfg) -> str:
    return _NOT_FOUND_LITERAL if cfg.write_not_found_literal else ""


def _empty_pubchem_record(cfg: PubChemCfg) -> dict[str, str]:
    missing = _missing_literal(cfg)
    return {
        "pubchem_cid": missing,
        "pubchem_iupac_name": missing,
        "pubchem_molecular_formula": missing,
        "pubchem_isomeric_smiles": missing,
        "pubchem_canonical_smiles": missing,
        "pubchem_inchi": missing,
        "pubchem_inchikey": missing,
    }


def _is_truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes"}
    if isinstance(value, (int, float)):
        return bool(value)
    try:
        return bool(value)
    except TypeError:  # pandas.NA raises TypeError
        return False


def _unique(values: list[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            result.append(value)
    return result


def _build_url(identifier: str, value: str, cfg: PubChemCfg) -> str | None:
    base = cfg.base.rstrip("/")
    builder = _IDENTIFIER_PATHS.get(identifier)
    if not builder:
        return None
    return builder(base, value)


def _request_with_backoff(
    url: str,
    cfg: PubChemCfg,
    *,
    deadline: float | None,
    now: Callable[[], float],
    sleeper: Callable[[float], None],
) -> dict[str, Any] | None:
    global _CACHE
    if _CACHE is None or _CACHE.ttl != cfg.cache_ttl:
        _CACHE = TTLCache(maxsize=1024, ttl=cfg.cache_ttl)

    cached = _CACHE.get(url)
    if cached is not None:
        return cast(dict[str, Any], cached)

    wait = cfg.backoff_initial_seconds or cfg.delay
    attempts = max(1, cfg.retries)

    for attempt in range(1, attempts + 1):
        if deadline is not None and now() >= deadline:
            return None

        get_limiter("pubchem", cfg.rps, cfg.burst).acquire()
        try:
            response = _session.get(
                url, timeout=(cfg.timeout_connect, cfg.timeout_read)
            )
        except requests.RequestException:
            status = None
        else:
            status = response.status_code
            if status == 404:
                return None
            if status == 429 or 500 <= status < 600:
                data = None
            elif status is not None and status >= 400:
                return None
            else:
                try:
                    data = cast(dict[str, Any], response.json())
                except ValueError:
                    return None
                _CACHE[url] = data
                return data

        if attempt >= attempts:
            break
        if deadline is not None and now() >= deadline:
            break
        if wait and wait > 0:
            sleeper(wait)
            wait *= 2

    return None


def _parse_properties(data: dict[str, Any], cfg: PubChemCfg) -> dict[str, str] | None:
    properties = data.get("PropertyTable", {}).get("Properties", [])
    if not properties:
        return None
    entry = properties[0]
    cid_raw = entry.get("CID")
    cid = validate_cid(str(cid_raw)) if cid_raw is not None else None
    if not cid:
        return None

    missing = _missing_literal(cfg)

    def field(name: str) -> str:
        value = entry.get(name)
        if value is None or value == "":
            return missing
        return str(value)

    return {
        "pubchem_cid": cid,
        "pubchem_iupac_name": field("IUPACName"),
        "pubchem_molecular_formula": field("MolecularFormula"),
        "pubchem_isomeric_smiles": field("IsomericSMILES"),
        "pubchem_canonical_smiles": field("CanonicalSMILES"),
        "pubchem_inchi": field("InChI"),
        "pubchem_inchikey": field("InChIKey"),
    }


def _resolve_with_identifier(
    identifier: str,
    value: str,
    cfg: PubChemCfg,
    *,
    deadline: float | None,
    now: Callable[[], float],
    sleeper: Callable[[float], None],
) -> dict[str, str] | None:
    url = _build_url(identifier, value, cfg)
    if not url:
        return None
    data = _request_with_backoff(url, cfg, deadline=deadline, now=now, sleeper=sleeper)
    if not data:
        return None
    return _parse_properties(data, cfg)


def _cid_candidates(record: Mapping[str, Any]) -> list[str]:
    raw = record.get("pubchem_cid")
    if not raw:
        return []
    parts = str(raw).split("|")
    return _unique([validate_cid(part.strip()) or "" for part in parts])


def _smiles_candidates(record: Mapping[str, Any], cfg: PubChemCfg) -> list[str]:
    candidates: list[str] = []

    def add(value: Any) -> None:
        if isinstance(value, str):
            text = value.strip()
            if text:
                candidates.append(text)

    if cfg.prefer_local_smiles:
        add(record.get("canonical_smiles"))
        add(record.get("smiles"))
        add(record.get("standard_smiles"))
        add(record.get("molecule_structures__canonical_smiles"))
    else:
        add(record.get("pubchem_canonical_smiles"))
        add(record.get("pubchem_isomeric_smiles"))
        add(record.get("canonical_smiles"))

    if cfg.use_parent_for_salts and record.get("salt_chembl_id"):
        add(record.get("parent_canonical_smiles"))
        add(record.get("parent_isomeric_smiles"))

    return _unique(candidates)


def _inchikey_candidates(record: Mapping[str, Any]) -> list[str]:
    value = record.get("standard_inchi_key")
    if isinstance(value, str):
        value = value.strip()
        if value:
            return [value]
    return []


def _inchi_candidates(record: Mapping[str, Any]) -> list[str]:
    value = record.get("standard_inchi")
    if isinstance(value, str):
        value = value.strip()
        if value:
            return [value]
    return []


def _name_candidates(record: Mapping[str, Any], cfg: PubChemCfg) -> list[str]:
    candidates: list[str] = []

    def add(value: Any) -> None:
        if isinstance(value, str):
            text = value.strip()
            if text:
                candidates.append(text)

    if cfg.use_parent_for_salts and record.get("salt_chembl_id"):
        add(record.get("parent_pref_name"))
        add(record.get("parent_name"))
        add(record.get("parent_molecule_pref_name"))

    add(record.get("pref_name"))
    add(record.get("synonyms"))

    return _unique(candidates)


def _populate_cache(
    cache: MutableMapping[tuple[str, str], dict[str, str]],
    identifier: str,
    value: str,
    record: dict[str, str],
) -> None:
    cache[(identifier, value)] = record

    cid = record.get("pubchem_cid", "")
    if cid:
        cache.setdefault(("cid", cid), record)

    canonical = record.get("pubchem_canonical_smiles", "")
    if canonical:
        cache.setdefault(("smiles", canonical), record)

    isomeric = record.get("pubchem_isomeric_smiles", "")
    if isomeric:
        cache.setdefault(("smiles", isomeric), record)

    inchikey = record.get("pubchem_inchikey", "")
    if inchikey:
        cache.setdefault(("inchikey", inchikey), record)

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


def resolve_pubchem_record(
    record: Mapping[str, Any],
    cfg: PubChemCfg,
    *,
    cache: MutableMapping[tuple[str, str], dict[str, str]] | None = None,
    now: Callable[[], float] | None = None,
    sleeper: Callable[[float], None] | None = None,
) -> dict[str, str]:
    """Resolve a PubChem record using multiple identifiers.

    The function iterates over identifiers defined in ``cfg.resolve_order`` and
    queries the PubChem API until a match is found. Responses are cached in
    ``cache`` by identifier and CID to avoid duplicate lookups.

    Parameters
    ----------
    record:
        Mapping containing ChEMBL test item fields such as
        ``canonical_smiles`` and ``standard_inchi_key``.
    cfg:
        PubChem configuration.
    cache:
        Optional cache mapping ``(identifier, value)`` tuples to resolved
        PubChem records. When provided, results are reused across calls.
    now, sleeper:
        Optional callables used for timeout calculations and sleeping between
        retries. They default to :func:`time.monotonic` and
        :func:`library.rate_limiter.sleep` respectively and exist primarily for
        testing.

    Returns
    -------
    dict[str, str]
        Dictionary containing PubChem enrichment fields. When no identifier can
        be resolved, the dictionary contains empty strings or ``"Not Found"``
        depending on ``cfg.write_not_found_literal``.
    """

    if not cfg.enable:
        return _empty_pubchem_record(cfg)

    polymer_flag = record.get("polymer_flag")
    if not cfg.allow_polymer and _is_truthy(polymer_flag):
        return _empty_pubchem_record(cfg)

    molecule_type = record.get("molecule_type")
    if (
        not cfg.allow_polymer
        and isinstance(molecule_type, str)
        and molecule_type.strip().lower() == "polymer"
    ):
        return _empty_pubchem_record(cfg)

    now_fn = now or time.monotonic
    sleep_fn = sleeper or sleep
    deadline = now_fn() + cfg.timeout_seconds if cfg.timeout_seconds > 0 else None

    candidates: dict[str, list[str]] = {
        "cid": _cid_candidates(record),
        "smiles": _smiles_candidates(record, cfg),
        "inchikey": _inchikey_candidates(record),
        "inchi": _inchi_candidates(record),
        "pref_name": _name_candidates(record, cfg),
    }

    for identifier in cfg.resolve_order:
        values = candidates.get(identifier, [])
        for value in values:
            key = (identifier, value)
            if cache is not None and key in cache:
                return cache[key]

            resolved = _resolve_with_identifier(
                identifier,
                value,
                cfg,
                deadline=deadline,
                now=now_fn,
                sleeper=sleep_fn,
            )
            if resolved:
                if cache is not None:
                    _populate_cache(cache, identifier, value, resolved)
                return resolved

    return _empty_pubchem_record(cfg)


def process_compound(compound_name: str, cfg: PubChemCfg) -> dict[str, str]:
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
    "resolve_pubchem_record",
    "process_compound",
    "Properties",
]
