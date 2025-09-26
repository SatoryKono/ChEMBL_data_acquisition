"""PubChem API client utilities.

This module provides functions to interact with the PubChem REST API.
The implementation is a Python translation of a PowerQuery script.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, MutableMapping
from dataclasses import dataclass
from time import monotonic
from typing import Any, Hashable, cast
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
    """Make an HTTP GET request and return parsed JSON."""

    global _CACHE
    if _CACHE is None or _CACHE.ttl != cfg.cache_ttl:
        _CACHE = TTLCache(maxsize=1024, ttl=cfg.cache_ttl)

    cached = _CACHE.get(url)
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
            logger.warning(
                "request_timeout", url=url, attempt=attempt, rps=cfg.rps
            )
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
            event_name = "request_rate_limited" if status == 429 else "request_server_error"
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
        assert _CACHE is not None
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


def _select_primary_cid(
    candidates: str | None,
    *,
    cache_key: str | None,
    identifier: str,
    value: str | None,
) -> str | None:
    """Return the first CID from ``candidates`` logging alternatives."""

    if not candidates:
        return None
    cid_list = [cid.strip() for cid in str(candidates).split("|") if cid.strip()]
    if not cid_list:
        return None
    primary = cid_list[0]
    if len(cid_list) > 1:
        logger.info(
            "pubchem_multiple_cid",
            chembl_id=cache_key,
            identifier=identifier,
            value=value,
            cid=primary,
            alternatives=cid_list[1:],
        )
    return primary


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


@dataclass(frozen=True)
class PubChemResolution:
    """Result of resolving a PubChem record."""

    cid: str | None
    source: str | None
    status: int | None = None


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


def resolve_pubchem_record(
    identifiers: Mapping[str, str | None],
    cfg: PubChemCfg,
    *,
    cid_cache: MutableMapping[str, str | None] | None = None,
    cache_key: str | None = None,
    resolution_cache: MutableMapping[Hashable, PubChemResolution] | None = None,
    resolution_key: Hashable | None = None,
) -> PubChemResolution:
    """Resolve a PubChem CID according to ``cfg.resolve_order``."""

    if resolution_cache is not None and resolution_key is not None:
        cached_resolution = resolution_cache.get(resolution_key)
        if cached_resolution is not None:
            return cached_resolution

    def _remember(resolution: PubChemResolution) -> PubChemResolution:
        if resolution_cache is not None and resolution_key is not None:
            resolution_cache[resolution_key] = resolution
        if resolution.cid and cid_cache is not None and cache_key:
            existing = cid_cache.get(cache_key)
            if existing != resolution.cid:
                cid_cache[cache_key] = resolution.cid
        return resolution

    def _resolve_cache() -> PubChemResolution | None:
        if cid_cache is not None and cache_key:
            cached_value = cid_cache.get(cache_key)
            if isinstance(cached_value, str) and cached_value:
                cid = _select_primary_cid(
                    cached_value,
                    cache_key=cache_key,
                    identifier="cache",
                    value=cached_value,
                )
                if cid:
                    if cached_value != cid:
                        cid_cache[cache_key] = cid
                    return PubChemResolution(cid=cid, source="cache")
        existing_cid = identifiers.get("pubchem_cid")
        if existing_cid:
            cid = _select_primary_cid(
                existing_cid,
                cache_key=cache_key,
                identifier="pubchem_cid",
                value=existing_cid,
            )
            if cid:
                return PubChemResolution(cid=cid, source="pubchem_cid")
        return None

    def _attempt_candidates(
        resolver: Callable[[str, PubChemCfg], str | None],
        candidates: tuple[tuple[str, str | None], ...],
    ) -> PubChemResolution | None:
        for identifier, value in candidates:
            if not value:
                continue
            resolved = resolver(value, cfg)
            cid = _select_primary_cid(
                resolved,
                cache_key=cache_key,
                identifier=identifier,
                value=value,
            )
            if cid:
                return PubChemResolution(cid=cid, source=identifier)
        return None

    def _resolve_smiles() -> PubChemResolution | None:
        candidates = (
            ("canonical_smiles", identifiers.get("canonical_smiles")),
            ("pubchem_canonical_smiles", identifiers.get("pubchem_canonical_smiles")),
            ("isomeric_smiles", identifiers.get("isomeric_smiles")),
            ("pubchem_isomeric_smiles", identifiers.get("pubchem_isomeric_smiles")),
        )
        return _attempt_candidates(get_cid_from_smiles, candidates)

    def _resolve_inchikey() -> PubChemResolution | None:
        candidates = (
            ("standard_inchi_key", identifiers.get("standard_inchi_key")),
            ("pubchem_inchikey", identifiers.get("pubchem_inchikey")),
        )
        return _attempt_candidates(get_cid_from_inchikey, candidates)

    def _resolve_inchi() -> PubChemResolution | None:
        candidates = (
            ("standard_inchi", identifiers.get("standard_inchi")),
            ("pubchem_inchi", identifiers.get("pubchem_inchi")),
        )
        return _attempt_candidates(get_cid_from_inchi, candidates)

    def _resolve_name() -> PubChemResolution | None:
        name_value = identifiers.get("pref_name")
        if not name_value:
            return None
        for identifier, resolver in (
            ("pref_name", get_cid),
            ("pref_name_partial", get_all_cid),
        ):
            resolved = resolver(name_value, cfg)
            cid = _select_primary_cid(
                resolved,
                cache_key=cache_key,
                identifier=identifier,
                value=name_value,
            )
            if cid:
                return PubChemResolution(cid=cid, source=identifier)
        return None

    handlers: dict[str, Callable[[], PubChemResolution | None]] = {
        "cache": _resolve_cache,
        "smiles": _resolve_smiles,
        "inchikey": _resolve_inchikey,
        "inchi": _resolve_inchi,
        "pref_name": _resolve_name,
        "name": _resolve_name,
    }

    for stage in cfg.resolve_order:
        handler = handlers.get(stage.lower())
        if handler is None:
            raise ValueError(f"Unknown PubChem resolve order entry: {stage!r}")
        resolution = handler()
        if resolution and resolution.cid:
            return _remember(resolution)

    final = PubChemResolution(cid=None, source=None)
    if resolution_cache is not None and resolution_key is not None:
        resolution_cache[resolution_key] = final
    return final


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
    "PubChemResolution",
]
