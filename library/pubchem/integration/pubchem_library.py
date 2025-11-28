"""PubChem API client utilities.

This module provides functions to interact with the PubChem REST API.
The implementation is a Python translation of a PowerQuery script.
"""

from __future__ import annotations

from collections.abc import Callable, Hashable, Mapping, MutableMapping
from dataclasses import dataclass

from ..clients.pubchem import (
    Properties,
    PubChemServiceUnavailable,
    get_all_cid,
    get_cid,
    get_cid_from_inchi,
    get_cid_from_inchikey,
    get_cid_from_smiles,
    get_properties,
    get_standard_name,
    init_session,
    make_request,
    url_encode,
    validate_cid,
)
from ..common.log import logger
from ..config import PubChemCfg


def select_primary_cid(
    candidates: str | None,
    *,
    identifier: str,
    value: str | None,
    context_id: str | None = None,
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
            chembl_id=context_id,
            identifier=identifier,
            value=value,
            cid=primary,
            alternatives=cid_list[1:],
        )
    return primary


@dataclass(frozen=True)
class PubChemResolution:
    """Result of resolving a PubChem record."""

    cid: str | None
    source: str | None
    status: int | None = None
    temporary_failure: bool = False


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
                cid = select_primary_cid(
                    cached_value,
                    identifier="cache",
                    value=cached_value,
                    context_id=cache_key,
                )
                if cid:
                    if cached_value != cid:
                        cid_cache[cache_key] = cid
                    return PubChemResolution(cid=cid, source="cache")
        existing_cid = identifiers.get("pubchem_cid")
        if existing_cid:
            cid = select_primary_cid(
                existing_cid,
                identifier="pubchem_cid",
                value=existing_cid,
                context_id=cache_key,
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
            cid = select_primary_cid(
                resolved,
                identifier=identifier,
                value=value,
                context_id=cache_key,
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
            cid = select_primary_cid(
                resolved,
                identifier=identifier,
                value=name_value,
                context_id=cache_key,
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
        try:
            resolution = handler()
        except PubChemServiceUnavailable as exc:
            log_fields: dict[str, object] = {"stage": stage}
            if exc.outcome:
                log_fields["outcome"] = exc.outcome
            if exc.status is not None:
                log_fields["status"] = exc.status
            for key, value in exc.details.items():
                if key == "status":
                    continue
                if key not in log_fields:
                    log_fields[key] = value
            cooldown_remaining = exc.details.get("cooldown_remaining")
            log_level = "INFO" if cooldown_remaining is not None else "WARNING"
            logger.log(log_level, "pubchem_unavailable", **log_fields)
            resolution = PubChemResolution(
                cid=None,
                source=None,
                status=exc.status,
                temporary_failure=True,
            )
            if resolution_cache is not None and resolution_key is not None:
                resolution_cache[resolution_key] = resolution
            return resolution
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
    props = (
        get_properties(cid, cfg)
        if cid
        else Properties(None, None, None, None, None, None)
    )
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
    "get_cid_from_smiles",
    "get_cid_from_inchi",
    "get_cid_from_inchikey",
    "get_cid",
    "get_all_cid",
    "get_standard_name",
    "get_properties",
    "resolve_pubchem_record",
    "process_compound",
    "Properties",
    "PubChemResolution",
]
