"""Target retrieval helpers for the ChEMBL API."""

from __future__ import annotations

from typing import Any, Iterable

import logging

import pandas as pd

from .chembl_client import _chunked, request_json
from .mapper_library import map_chembl_to_uniprot
from .config import Config

logger = logging.getLogger(__name__)

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
]

EMPTY_TARGET: dict[str, str] = {field: "" for field in TARGET_FIELDS}


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


def _parse_uniprot_id(xrefs: list[dict[str, str]], chembl_id: str) -> tuple[str, str]:
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
        mapping_uniprot_id = map_chembl_to_uniprot(chembl_id) or ""
    except Exception as exc:  # pragma: no cover - network failure paths
        logger.warning("UniProt mapping request failed for %s: %s", chembl_id, exc)
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


def _parse_target_record(data: dict[str, Any]) -> dict[str, Any]:
    """Transform a raw target record into a flat dictionary."""
    components = _get_items(data.get("target_components"), "target_component")
    if not components:
        logger.debug("No components found in target record: %s", data)
        return dict(EMPTY_TARGET)

    comp = components[0]
    synonyms = _get_items(
        comp.get("target_component_synonyms"), "target_component_synonym"
    )
    xrefs = _get_items(comp.get("target_component_xrefs"), "target")

    gene_syn = _parse_gene_synonyms(synonyms)
    ec_code = _parse_ec_codes(synonyms)
    alt_name = _parse_alt_names(synonyms)
    uniprot_id, mapping_uniprot_id = _parse_uniprot_id(
        xrefs, data.get("target_chembl_id", "")
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
        }
    )
    return res


def get_target(
    chembl_target_id: str,
    *,
    cfg: Config,
    timeout: float | None = None,
) -> dict[str, str]:
    """Fetch target information for a single ChEMBL identifier."""
    if chembl_target_id in {"", "#N/A"}:
        return dict(EMPTY_TARGET)
    url = "https://www.ebi.ac.uk/chembl/api/data/target/{id}.json".format(
        id=chembl_target_id
    )
    data = request_json(cfg, url, timeout=timeout)
    target_list = _get_items(data, "target")
    if not target_list:
        return dict(EMPTY_TARGET)
    return _parse_target_record(target_list[0])


def get_targets(
    ids: Iterable[str],
    *,
    cfg: Config,
    chunk_size: int = 5,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Fetch target records for ``ids``."""
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=TARGET_FIELDS)

    records: list[dict[str, Any]] = []
    base = "https://www.ebi.ac.uk/chembl/api/data/target.json?format=json"
    for chunk in _chunked(valid, chunk_size):
        url = f"{base}&target_chembl_id__in={','.join(chunk)}"
        data = request_json(cfg, url, timeout=timeout)
        items = data.get("targets") or data.get("target") or []
        records.extend(_parse_target_record(item) for item in items)
    if not records:
        return pd.DataFrame(columns=TARGET_FIELDS)
    df = pd.DataFrame(records)
    return df.reindex(columns=TARGET_FIELDS)


def extend_target(
    df: pd.DataFrame, *, cfg: Config, id_column: str = "target_chembl_id"
) -> pd.DataFrame:
    """Augment ``df`` with columns returned from :func:`get_target`.

    Parameters
    ----------
    df:
        DataFrame containing a column with ChEMBL target identifiers.
    id_column:
        Name of the column holding the identifiers.

    """
    if id_column not in df.columns:
        raise ValueError(f"missing required column: {id_column}")
    targets = [get_target(i, cfg=cfg) for i in df[id_column].fillna("")]
    extra = pd.DataFrame(targets)
    return pd.concat([df.reset_index(drop=True), extra], axis=1)


__all__ = [
    "get_target",
    "get_targets",
    "extend_target",
    "TARGET_FIELDS",
]
