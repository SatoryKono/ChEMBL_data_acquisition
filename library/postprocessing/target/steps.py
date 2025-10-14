"""Reproducible target postprocessing pipeline steps."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import pandas as pd
import requests

from library.clients.chembl_client import ChemblClient
from library.clients.gtopdb_client import GtoPdbClient
from library.clients.uniprot_client import UniProtClient
from library.schemas.target_schema import TargetSchema
from library.utils.logging import get_logger
from library.utils.qc_report import TableQualityProfiler, build_reports_from_profiler

logger = get_logger(__name__)

_TARGET_COLUMNS: Sequence[str] = (
    "target_chembl_id",
    "pref_name",
    "target_type",
    "organism",
    "uniprot_id",
    "gene_symbol",
    "target_class",
    "protein_family",
    "synonyms",
)


@dataclass(slots=True)
class TargetData:
    """Container grouping the primary dataset and QC artefacts."""

    dataset: pd.DataFrame
    quality_report: pd.DataFrame
    correlation_report: pd.DataFrame
    qc_summary: dict[str, Any]


def fetch_normalize_target(
    limit: int,
    *,
    chembl_client: ChemblClient | None = None,
    uniprot_client: UniProtClient | None = None,
    gtopdb_client: GtoPdbClient | None = None,
) -> pd.DataFrame:
    """Return the validated target dataset enriched with external metadata."""

    if limit <= 0:
        raise ValueError("limit must be positive")

    chembl = chembl_client if chembl_client is not None else ChemblClient()
    logger.info("target_fetch_start", limit=limit)
    target_payload = _fetch_targets(chembl, limit)
    base_frame = _normalise_targets(target_payload)
    logger.info("target_fetch_success", rows=len(base_frame))

    if base_frame.empty:
        validated = TargetSchema.validate(base_frame)
        return validated

    uniprot_ids = _extract_uniprot_ids(base_frame)
    uni_client = uniprot_client if uniprot_client is not None else UniProtClient()
    gto_client = gtopdb_client if gtopdb_client is not None else GtoPdbClient()

    uniprot_metadata = _load_uniprot_metadata(uni_client, uniprot_ids)
    gtopdb_metadata = _load_gtopdb_metadata(gto_client, uniprot_ids)

    enriched = enrich_target_metadata(base_frame, uniprot_metadata, gtopdb_metadata)
    validated = TargetSchema.validate(enriched)
    return validated


def enrich_target_metadata(
    frame: pd.DataFrame,
    uniprot_metadata: pd.DataFrame,
    gtopdb_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Merge UniProt and GtoPdb annotations without overwriting curated data."""

    working = frame.copy()
    for column in ("target_class", "protein_family", "synonyms"):
        if column not in working.columns:
            working[column] = pd.Series([pd.NA] * len(working), dtype="string")
    working["uniprot_id"] = working["uniprot_id"].astype("string")

    merged = _combine_metadata(working, gtopdb_metadata, keys=("target_class",))
    merged = _combine_metadata(merged, uniprot_metadata, keys=("protein_family", "synonyms"))
    return merged


def generate_target_reports(frame: pd.DataFrame) -> TargetData:
    """Compute QC and correlation artefacts for ``frame``."""

    profiler = TableQualityProfiler()
    profiler.consume(frame)
    quality_report, correlation_report = build_reports_from_profiler(profiler)
    qc_summary = {
        "row_count": int(frame.shape[0]),
        "column_count": int(frame.shape[1]),
        "non_null_ratio": float(frame.notna().mean().mean()) if not frame.empty else 0.0,
    }
    return TargetData(
        dataset=frame,
        quality_report=quality_report,
        correlation_report=correlation_report,
        qc_summary=qc_summary,
    )


def _fetch_targets(client: ChemblClient, limit: int) -> list[Mapping[str, Any]]:
    records: list[Mapping[str, Any]] = []
    offset = 0
    while len(records) < limit:
        page_limit = min(200, limit - len(records))
        payload = client.list_targets(limit=page_limit, offset=offset)
        items = _extract_payload_targets(payload)
        if not items:
            break
        records.extend(items)
        offset += len(items)
        page_meta = payload.get("page_meta") if isinstance(payload, Mapping) else None
        next_url = page_meta.get("next") if isinstance(page_meta, Mapping) else None
        if not next_url:
            break
    return records[:limit]


def _extract_payload_targets(payload: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    for key in ("targets", "items", "data", "results"):
        value = payload.get(key)
        if isinstance(value, list):
            return [item for item in value if isinstance(item, Mapping)]
    return []


def _normalise_targets(records: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for entry in records:
        target_id = str(entry.get("target_chembl_id") or "").strip()
        if not target_id:
            continue
        pref_name = _clean_string(entry.get("pref_name"))
        target_type = _clean_string(entry.get("target_type"))
        organism = _clean_string(entry.get("organism"))
        components = entry.get("target_components")
        if not isinstance(components, Sequence):
            components = []
        accession = _first_non_empty(
            [entry.get("accession"), *_extract_accessions(components)]
        )
        if accession is None:
            logger.warning("target_missing_uniprot", target_chembl_id=target_id)
            continue
        gene_symbol = _derive_gene_symbol(entry, components)
        synonyms = _collect_synonyms(entry, components)
        rows.append(
            {
                "target_chembl_id": target_id,
                "pref_name": pref_name or target_id,
                "target_type": target_type or "unknown",
                "organism": organism or "unknown",
                "uniprot_id": str(accession),
                "gene_symbol": gene_symbol or target_id,
                "target_class": _clean_string(entry.get("target_class")),
                "protein_family": pd.NA,
                "synonyms": synonyms,
            }
        )
    frame = pd.DataFrame(rows, columns=_TARGET_COLUMNS)
    if frame.empty:
        return frame
    for column in _TARGET_COLUMNS:
        frame[column] = frame[column].astype("string")
    frame = frame.drop_duplicates(subset=["target_chembl_id"]).reset_index(drop=True)
    return frame


def _extract_uniprot_ids(frame: pd.DataFrame) -> list[str]:
    ids = frame["uniprot_id"].dropna().astype("string")
    unique_ids = sorted({str(value) for value in ids if value != "<NA>"})
    return unique_ids


def _load_uniprot_metadata(client: UniProtClient, uniprot_ids: Sequence[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for accession in uniprot_ids:
        try:
            payload = client.get_entry(accession)
        except requests.exceptions.RequestException as exc:
            logger.warning("uniprot_fetch_failed", uniprot_id=accession, error=str(exc))
            continue
        family = _extract_uniprot_family(payload)
        synonyms = _normalise_synonym_list(_extract_uniprot_synonyms(payload))
        rows.append(
            {
                "uniprot_id": str(accession),
                "protein_family": family,
                "synonyms": synonyms,
            }
        )
    frame = pd.DataFrame(rows, columns=["uniprot_id", "protein_family", "synonyms"])
    if frame.empty:
        return frame
    frame["uniprot_id"] = frame["uniprot_id"].astype("string")
    frame["protein_family"] = frame["protein_family"].astype("string")
    frame["synonyms"] = frame["synonyms"].astype("string")
    return frame


def _load_gtopdb_metadata(client: GtoPdbClient, uniprot_ids: Sequence[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for accession in uniprot_ids:
        try:
            payload = client.get_target_by_uniprot(accession)
        except requests.exceptions.RequestException as exc:
            logger.warning("gtopdb_fetch_failed", uniprot_id=accession, error=str(exc))
            continue
        target_class = _extract_gtopdb_class(payload)
        rows.append({"uniprot_id": str(accession), "target_class": target_class})
    frame = pd.DataFrame(rows, columns=["uniprot_id", "target_class"])
    if frame.empty:
        return frame
    frame["uniprot_id"] = frame["uniprot_id"].astype("string")
    frame["target_class"] = frame["target_class"].astype("string")
    return frame


def _combine_metadata(
    frame: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    keys: Sequence[str],
) -> pd.DataFrame:
    if metadata.empty:
        return frame
    working = frame.copy()
    meta = metadata.copy()
    meta["uniprot_id"] = meta["uniprot_id"].astype("string")
    combined = working.merge(meta, on="uniprot_id", how="left", suffixes=("", "_meta"))
    for key in keys:
        meta_column = f"{key}_meta"
        if meta_column in combined.columns:
            combined[key] = combined[key].combine_first(combined[meta_column])
            combined = combined.drop(columns=[meta_column])
    return combined


def _extract_accessions(components: Sequence[Mapping[str, Any]]) -> list[str]:
    values: list[str] = []
    for component in components:
        accession = component.get("accession") if isinstance(component, Mapping) else None
        if isinstance(accession, str) and accession.strip():
            values.append(accession.strip())
    return values


def _derive_gene_symbol(
    entry: Mapping[str, Any], components: Sequence[Mapping[str, Any]]
) -> str | None:
    candidates: list[str | None] = []
    if isinstance(entry.get("gene_symbol"), str):
        candidates.append(entry.get("gene_symbol"))
    for component in components:
        if not isinstance(component, Mapping):
            continue
        for key in ("gene_symbol", "component_gene_symbol", "gene_name"):
            value = component.get(key)
            if isinstance(value, str):
                candidates.append(value)
        synonyms = component.get("component_synonyms")
        if isinstance(synonyms, Sequence):
            for synonym in synonyms:
                if not isinstance(synonym, Mapping):
                    continue
                syn_type = synonym.get("syn_type") or synonym.get("type")
                if isinstance(syn_type, str) and syn_type.upper() in {"GENE_SYMBOL", "GENE"}:
                    syn_value = synonym.get("synonym")
                    if isinstance(syn_value, str):
                        candidates.append(syn_value)
    for candidate in candidates:
        cleaned = _clean_string(candidate)
        if cleaned:
            return cleaned
    return None


def _collect_synonyms(
    entry: Mapping[str, Any], components: Sequence[Mapping[str, Any]]
) -> str | None:
    values: list[str] = []
    entry_synonyms = entry.get("synonyms")
    if isinstance(entry_synonyms, Sequence):
        for synonym in entry_synonyms:
            cleaned = _clean_string(synonym)
            if cleaned:
                values.append(cleaned)
    for component in components:
        if not isinstance(component, Mapping):
            continue
        synonyms = component.get("component_synonyms")
        if isinstance(synonyms, Sequence):
            for synonym in synonyms:
                if isinstance(synonym, Mapping):
                    cleaned = _clean_string(synonym.get("synonym"))
                else:
                    cleaned = _clean_string(synonym)
                if cleaned:
                    values.append(cleaned)
    normalised = _normalise_synonym_list(values)
    return normalised


def _extract_uniprot_family(payload: Any) -> str | None:
    candidates = []
    if isinstance(payload, Mapping):
        for key in ("protein_families", "families", "proteinFamily", "family"):
            value = payload.get(key)
            candidates.extend(_unwrap_strings(value))
        keywords = payload.get("keywords")
        candidates.extend(_unwrap_strings(keywords))
    if candidates:
        return _clean_string(candidates[0])
    return None


def _extract_uniprot_synonyms(payload: Any) -> list[str]:
    synonyms: list[str] = []
    if isinstance(payload, Mapping):
        for key in ("synonyms", "alternative_names", "alternativeNames", "gene_names"):
            synonyms.extend(_unwrap_strings(payload.get(key)))
        gene_section = payload.get("genes")
        if isinstance(gene_section, Sequence):
            for gene in gene_section:
                if isinstance(gene, Mapping):
                    synonyms.extend(_unwrap_strings(gene.get("synonyms")))
                    synonyms.extend(_unwrap_strings(gene.get("orfNames")))
    return [value for value in (_clean_string(item) for item in synonyms) if value]


def _extract_gtopdb_class(payload: Any) -> str | None:
    if isinstance(payload, Mapping):
        for key in ("target_class", "targetClass", "type", "class"):
            candidate = _clean_string(payload.get(key))
            if candidate:
                return candidate
        entries = payload.get("targets")
        if isinstance(entries, Sequence):
            for entry in entries:
                candidate = _extract_gtopdb_class(entry)
                if candidate:
                    return candidate
    if isinstance(payload, Sequence):
        for item in payload:
            candidate = _extract_gtopdb_class(item)
            if candidate:
                return candidate
    return None


def _normalise_synonym_list(values: Iterable[str]) -> str | None:
    unique = list(dict.fromkeys(_clean_string(value) for value in values if value))
    if not unique:
        return None
    return "; ".join(unique)


def _unwrap_strings(value: Any) -> list[str]:
    if isinstance(value, str):
        return [value]
    if isinstance(value, Mapping):
        collected: list[str] = []
        for maybe in value.values():
            collected.extend(_unwrap_strings(maybe))
        return collected
    if isinstance(value, Sequence):
        collected = []
        for item in value:
            collected.extend(_unwrap_strings(item))
        return collected
    return []


def _clean_string(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    cleaned = value.strip()
    return cleaned or None


def _first_non_empty(values: Iterable[Any]) -> str | None:
    for value in values:
        cleaned = _clean_string(value)
        if cleaned:
            return cleaned
    return None


__all__ = [
    "TargetData",
    "fetch_normalize_target",
    "generate_target_reports",
    "enrich_target_metadata",
]
