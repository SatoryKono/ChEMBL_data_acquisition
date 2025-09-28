"""Target retrieval helpers for the ChEMBL API."""

from __future__ import annotations

import json
import re
from collections.abc import Iterable
from typing import Any

import pandas as pd

from ..clients import ChemblClient, chunked
from ..clients.mapper import map_chembl_to_uniprot
from ..config import ApiCfg, UniprotMappingCfg
from ..utils.logging import logger

TARGET_FIELDS = [
    "pref_name",
    "target_chembl_id",
    "component_description",
    "component_id",
    "relationship",
    "gene",
    "uniprot_id",
    "mapping_uniprot_id",
    "chembl_alternative_name",
    "ec_code",
    "hgnc_name",
    "hgnc_id",
    "target_type",
    "tax_id",
    "species_group_flag",
    "target_components",
    "protein_classifications",
    "cross_references",
    "reaction_ec_numbers",
]

EMPTY_TARGET: dict[str, str] = {field: "" for field in TARGET_FIELDS}

TARGET_INCLUDE_PARAMS = "protein_classifications,cross_references"

REACTION_EC_EXCLUDED_XREFS = {
    "REACTOME",
    "RHEA",
    "METACYC",
    "EC_REACTION",
}

_EC_TOKEN_SPLIT = re.compile(r"[|;,/\\\s]+")
_EC_FULL_PATTERN = re.compile(r"^\d+(?:\.(?:\d+|-)){3}$")


def normalize_reaction_ec_numbers(values: Iterable[str | None]) -> str:
    """Return a pipe-delimited string of sanitized EC numbers from ``values``."""

    numbers = _collect_normalized_ec_tokens(values)
    return "|".join(sorted(numbers))


def _collect_normalized_ec_tokens(values: Iterable[str | None]) -> set[str]:
    """Return a set of EC numbers extracted from ``values``."""

    numbers: set[str] = set()
    for value in values:
        if not isinstance(value, str) or not value:
            continue
        for token in _split_ec_candidates(value):
            cleaned = _normalise_ec_token(token)
            if cleaned and _EC_FULL_PATTERN.fullmatch(cleaned):
                numbers.add(cleaned)
    return numbers


def _split_ec_candidates(value: str) -> list[str]:
    """Return tokens parsed from ``value`` using standard separators."""

    return [token for token in _EC_TOKEN_SPLIT.split(value) if token]


def _normalise_ec_token(token: str) -> str:
    """Return a cleaned EC token stripped of prefixes and whitespace."""

    token = token.strip()
    if not token:
        return ""
    upper = token.upper()
    if upper.startswith("EC"):
        token = token[2:]
        token = token.lstrip(":._- ")
    return token.strip()


def _parse_gene_synonyms(synonyms: list[dict[str, str]]) -> str:
    """Return a sorted, pipe separated list of gene synonyms."""
    names = {
        s["component_synonym"]
        for s in synonyms
        if s.get("syn_type") in {"GENE_SYMBOL", "GENE_SYMBOL_OTHER"}
    }
    return "|".join(sorted(names))


def _parse_ec_codes(synonyms: list[dict[str, str]]) -> str:
    """Return a sorted, pipe separated list of EC numbers."""
    codes = {
        s["component_synonym"] for s in synonyms if s.get("syn_type") == "EC_NUMBER"
    }
    return "|".join(sorted(codes))


def _parse_alt_names(synonyms: list[dict[str, str]]) -> str:
    """Return a sorted, pipe separated list of UniProt alternative names."""
    names = {s["component_synonym"] for s in synonyms if s.get("syn_type") == "UNIPROT"}
    return "|".join(sorted(names))


def _parse_uniprot_id(
    xrefs: list[dict[str, str]], chembl_id: str, mapping_cfg: UniprotMappingCfg
) -> tuple[str, str]:
    """Return UniProt IDs from cross references and mapping."""
    uniprot_id = ""
    for x in xrefs:
        src = (x.get("xref_src_db") or "").upper()
        if src in {"UNIPROT", "UNIPROT ACCESSION", "UNIPROT ACC", "UNIPROTKB"}:
            ident = x.get("xref_id", "")
            if ident:
                uniprot_id = ident
                break
    try:
        mapping_uniprot_id = map_chembl_to_uniprot(chembl_id, mapping_cfg) or ""
    except Exception as exc:  # pragma: no cover - network failure paths
        logger.warning(
            "uniprot_mapping_error",
            chembl_id=str(chembl_id),
            error=str(exc),
        )
        mapping_uniprot_id = ""
    return uniprot_id, mapping_uniprot_id


def _parse_hgnc(xrefs: list[dict[str, str]]) -> tuple[str, str]:
    """Extract HGNC name and identifier from a list of cross references."""
    for x in xrefs:
        if x.get("xref_src_db") == "HGNC":
            name = x.get("xref_name", "")
            ident = x.get("xref_id", "")
            hgnc_id = ident.split(":")[-1] if ident else ""
            return name, hgnc_id
    return "", ""


def _stringify(value: Any) -> str:
    """Return a string representation of ``value`` avoiding ``None``."""

    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _serialize_structure(value: Any) -> str:
    """Return a JSON string for ``value`` or ``""`` when empty."""

    if value in (None, "", [], {}):
        return ""
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    except (TypeError, ValueError):
        return _stringify(value)


def _collect_reaction_ec_numbers(components: list[dict[str, Any]]) -> str:
    """Return pipe-delimited reaction EC numbers discovered in ``components``."""

    candidates: list[str] = []
    for component in components:
        synonyms = _get_items(
            component.get("target_component_synonyms"), "target_component_synonym"
        )
        for synonym in synonyms:
            syn_type = (synonym.get("syn_type") or "").upper()
            if syn_type in {"REACTION", "REACTION_NUMBER", "EC_REACTION_NUMBER"}:
                value = synonym.get("component_synonym", "")
                if isinstance(value, str) and value:
                    candidates.append(value)
        xrefs = _get_items(component.get("target_component_xrefs"), "target")
        for xref in xrefs:
            src = (xref.get("xref_src_db") or "").upper()
            if src in REACTION_EC_EXCLUDED_XREFS:
                continue
            value = xref.get("xref_id")
            if isinstance(value, str) and value:
                candidates.append(value)
    return normalize_reaction_ec_numbers(candidates)


def _get_items(container: Any, key: str) -> list[Any]:
    """Return a list of items from a container that may be a dict or list."""
    if isinstance(container, dict):
        items = container.get(key, [])
    else:
        items = container or []
    if isinstance(items, dict):
        return [items]
    if isinstance(items, list):
        return items
    return []


def _extract_target_payload(data: Any) -> dict[str, Any]:
    """Return a target record from ``data`` regardless of envelope shape."""

    if isinstance(data, dict):
        if "target" in data:
            items = data["target"]
            if isinstance(items, list):
                return items[0] if items else {}
            if isinstance(items, dict):
                return items
        if "targets" in data:
            items = data["targets"]
            if isinstance(items, list):
                return items[0] if items else {}
            if isinstance(items, dict):
                return items
        return data
    if isinstance(data, list):
        return data[0] if data else {}
    return {}


def _parse_target_record(
    data: dict[str, Any], mapping_cfg: UniprotMappingCfg
) -> dict[str, Any]:
    """Transform a raw target record into a flat dictionary."""
    components = _get_items(data.get("target_components"), "target_component")
    if not components:
        logger.debug("No components found in target record: %s", data)
        components = []
    comp = components[0] if components else {}
    synonyms = _get_items(
        comp.get("target_component_synonyms"), "target_component_synonym"
    )
    xrefs = _get_items(comp.get("target_component_xrefs"), "target")

    gene_syn = _parse_gene_synonyms(synonyms)
    ec_code = _parse_ec_codes(synonyms)
    alt_name = _parse_alt_names(synonyms)
    uniprot_id, mapping_uniprot_id = _parse_uniprot_id(
        xrefs, data.get("target_chembl_id", ""), mapping_cfg
    )
    hgnc_name, hgnc_id = _parse_hgnc(xrefs)

    res = dict(EMPTY_TARGET)
    res.update(
        {
            "pref_name": data.get("pref_name", ""),
            "target_chembl_id": data.get("target_chembl_id", ""),
            "component_description": comp.get("component_description", ""),
            "component_id": str(comp.get("component_id", "")),
            "relationship": data.get("target_type", ""),
            "gene": gene_syn,
            "uniprot_id": uniprot_id,
            "mapping_uniprot_id": mapping_uniprot_id,
            "chembl_alternative_name": alt_name,
            "ec_code": ec_code,
            "hgnc_name": hgnc_name,
            "hgnc_id": hgnc_id,
            "target_type": _stringify(data.get("target_type")),
            "tax_id": _stringify(data.get("tax_id")),
            "species_group_flag": _stringify(data.get("species_group_flag")),
            "target_components": _serialize_structure(components),
            "protein_classifications": _serialize_structure(
                data.get("protein_classifications")
            ),
            "cross_references": _serialize_structure(data.get("cross_references")),
            "reaction_ec_numbers": _collect_reaction_ec_numbers(components),
        }
    )
    return res


def get_target(
    chembl_target_id: str,
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    mapping_cfg: UniprotMappingCfg,
    timeout: float | None = None,
) -> dict[str, str]:
    """Fetch target information for a single ChEMBL identifier.

    Parameters
    ----------
    chembl_target_id:
        ChEMBL target identifier (e.g., ``"CHEMBL203"``).
    cfg:
        API configuration providing base URL and timeouts.
    client:
        ChemblClient used for HTTP requests and caching.
    mapping_cfg:
        Configuration for the UniProt mapping service.
    timeout:
        Optional override for the HTTP request timeout.
    """
    if chembl_target_id in {"", "#N/A"}:
        return dict(EMPTY_TARGET)
    base = cfg.chembl_base.rstrip("/")
    url = f"{base}/target/{chembl_target_id}.json?include={TARGET_INCLUDE_PARAMS}"
    effective_timeout = timeout if timeout is not None else cfg.timeout_read
    data = client.request_json(url, cfg=cfg, timeout=effective_timeout)
    payload = _extract_target_payload(data)
    if not payload:
        return dict(EMPTY_TARGET)
    return _parse_target_record(payload, mapping_cfg)


def get_targets(
    ids: Iterable[str],
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    mapping_cfg: UniprotMappingCfg,
    chunk_size: int = 5,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Fetch target records for ``ids``.

    Parameters
    ----------
    ids:
        ChEMBL target identifiers to retrieve.
    cfg:
        API configuration with base URL and timeouts.
    client:
        ChemblClient used for HTTP requests and caching.
    mapping_cfg:
        Settings for the UniProt mapping service.
    chunk_size:
        Number of identifiers to request per HTTP call.
    timeout:
        Optional override for the HTTP request timeout.
    """
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=TARGET_FIELDS)

    records: list[dict[str, Any]] = []
    base = (
        f"{cfg.chembl_base.rstrip('/')}/target.json?format=json"
        f"&include={TARGET_INCLUDE_PARAMS}"
    )
    effective_timeout = timeout if timeout is not None else cfg.timeout_read
    for chunk in chunked(valid, chunk_size):
        url = f"{base}&target_chembl_id__in={','.join(chunk)}"
        data = client.request_json(url, cfg=cfg, timeout=effective_timeout)
        items = data.get("targets") or data.get("target") or []
        records.extend(_parse_target_record(item, mapping_cfg) for item in items)
    if not records:
        return pd.DataFrame(columns=TARGET_FIELDS)
    df = pd.DataFrame(records)
    return df.reindex(columns=TARGET_FIELDS)


def extend_target(
    df: pd.DataFrame,
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    mapping_cfg: UniprotMappingCfg,
    id_column: str = "target_chembl_id",
) -> pd.DataFrame:
    """Augment ``df`` with columns returned from :func:`get_target`.

    Parameters
    ----------
    df:
        DataFrame containing a column with ChEMBL target identifiers.
    cfg:
        API configuration providing base URL and timeouts.
    client:
        ChemblClient used for HTTP requests and caching.
    mapping_cfg:
        Settings for the UniProt mapping service.
    id_column:
        Name of the column holding the identifiers.

    """
    if id_column not in df.columns:
        raise ValueError(f"missing required column: {id_column}")
    targets = [
        get_target(i, cfg=cfg, client=client, mapping_cfg=mapping_cfg)
        for i in df[id_column].fillna("")
    ]
    extra = pd.DataFrame(targets)
    return pd.concat([df.reset_index(drop=True), extra], axis=1)


__all__ = [
    "get_target",
    "get_targets",
    "extend_target",
    "TARGET_FIELDS",
]
