"""Utility functions for extracting information from UniProt data.

This module consolidates the logic previously spread across
``uniprot_names.py`` and ``uniprot_batch_names.py`` into a reusable
library.  It provides helpers to parse UniProt JSON structures, gather
protein and gene names, organism taxonomy, and batch-process multiple
entries from a CSV file.

The most commonly used functions are:

``fetch_uniprot(uniprot_id)``
    Retrieve a UniProt JSON entry from the REST API given an accession ID.

``extract_names(data)``
    Parse a UniProt JSON object and return a set of all protein and gene
    names found in the entry.

``extract_organism(data)``
    Extract genus, superkingdom, phylum and taxon ID information from a
    UniProt JSON object.  Returns a dictionary with these fields.

``iter_ids(csv_path)``
    Read a CSV file containing a ``uniprot_id`` column and yield each ID.

``collect_info(uid, data_dir="uniprot")``
    Given a UniProt accession and directory containing ``<uid>.json`` files
    (default: ``"uniprot"``), return a dictionary with the accession, all
    names, and organism taxonomy data.

``process(input_csv, output_csv, data_dir="uniprot")``
    Batch-process a CSV of UniProt IDs and write an output CSV with names and
    organism information for each ID.
"""

from __future__ import annotations

import csv
import json
import time
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import requests

from ..clients.uniprot import (
    UniProtFetchError,
    fetch_uniprot,
)
from ..clients.uniprot import (
    get_session as get_uniprot_session,
)
from ..clients.uniprot import (
    init_session as init_uniprot_session,
)
from ..common.log import logger
from ..common.rate_limiter import get_limiter
from ..config import IupharCfg, UniprotCfg

_DEFAULT_UNIPROT_DATA_DIR = Path("uniprot")


init_session = init_uniprot_session


__all__ = [
    "init_session",
    "fetch_uniprot",
    "extract_names",
    "extract_uniprotkb_id",
    "extract_secondary_accessions",
    "extract_names_for_secondary_accessions",
    "extract_recommended_name",
    "extract_gene_name",
    "extract_keywords",
    "extract_isoform",
    "extract_crossrefs",
    "extract_ptm",
    "extract_activity",
    "extract_organism",
    "iter_ids",
    "collect_info",
    "process",
    "UniProtFetchError",
]


def _canonical_accession(accession: str) -> str | None:
    """Return canonical accession for ``accession`` when it targets an isoform."""

    if not isinstance(accession, str):
        return None
    head, sep, tail = accession.partition("-")
    if not sep:
        return None
    head = head.strip()
    tail = tail.strip()
    if not head or not tail:
        return None
    return head


UNIPROT_OUTPUT_COLUMNS: list[str] = [
    "uniprot_id",
    "names",
    "genus",
    "superkingdom",
    "phylum",
    "lineage_class",
    "taxon_id",
    "sequence_length",
    "molecular_function",
    "cellular_component",
    "ec_numbers",
    "subcellular_location",
    "topology",
    "transmembrane",
    "intramembrane",
    "glycosylation",
    "lipidation",
    "disulfide_bond",
    "modified_residue",
    "phosphorylation",
    "acetylation",
    "ubiquitination",
    "signal_peptide",
    "propeptide",
    "isoform_names",
    "isoform_ids",
    "isoform_synonyms",
    "GuidetoPHARMACOLOGY",
    "family",
    "SUPFAM",
    "PROSITE",
    "InterPro",
    "Pfam",
    "PRINTS",
    "TCDB",
    "reactions",
    "reaction_ec_numbers",
    "gtop_natural_ligands_n",
    "gtop_interactions_n",
    "gtop_function_text_short",
    "xref_pdb",
    "xref_alphafold",
    "xref_ensembl",
    "uniprot_last_update",
    "uniprot_version",
    "pipeline_version",
    "timestamp_utc",
    "uniProtkbId",
    "uniProtkbIdFallback",
    "secondaryAccessions",
    "recommendedName",
    "geneName",
    "secondaryAccessionNames",
]


_GTOP_JSON_FAILURE_CACHE: set[tuple[str, str]] = set()
_GTOP_NON_JSON_CONTENT_TYPE_CACHE: set[tuple[str, str]] = set()
_GTOP_SKIPPED_FAILURE_LOG: set[tuple[str, str]] = set()


@dataclass
class _CircuitState:
    """State container tracking failures per Guide-to-Pharmacology endpoint."""

    failure_count: int = 0
    open: bool = False
    opened_until: float = 0.0
    last_holdoff: float = 0.0
    reset_pending: bool = False


@dataclass
class _CircuitDecision:
    """Outcome returned when consulting the circuit breaker before a request."""

    allow_call: bool
    remaining: float = 0.0
    reset_holdoff: float | None = None


@dataclass
class _CircuitBreaker:
    """Circuit breaker that throttles failing Guide-to-Pharmacology endpoints."""

    _states: dict[str, _CircuitState] = field(default_factory=dict)
    last_opened_at: float = 0.0
    last_holdoff: float = 0.0
    last_opened_endpoint: str | None = None

    def reset(self) -> None:
        """Reset the breaker state."""

        self._states.clear()
        self.last_opened_at = 0.0
        self.last_holdoff = 0.0
        self.last_opened_endpoint = None

    def before_call(self, endpoint: str) -> _CircuitDecision:
        """Return a decision describing whether a request should be attempted."""

        state = self._states.get(endpoint)
        if state is None:
            return _CircuitDecision(True)

        now = time.monotonic()
        if state.open:
            if now < state.opened_until:
                remaining = state.opened_until - now
                return _CircuitDecision(
                    False, remaining=remaining if remaining > 0 else 0.0
                )
            state.open = False
            state.reset_pending = True

        if state.reset_pending:
            state.reset_pending = False
            return _CircuitDecision(True, reset_holdoff=state.last_holdoff)

        return _CircuitDecision(True)

    def record_success(self, endpoint: str) -> None:
        """Clear any accumulated failures for ``endpoint`` after success."""

        if endpoint in self._states:
            del self._states[endpoint]

    def record_failure(self, endpoint: str, holdoff: float) -> bool:
        """Record a failure for ``endpoint`` and open the breaker if needed."""

        state = self._states.setdefault(endpoint, _CircuitState())
        state.failure_count += 1
        if state.failure_count < _GTOP_CIRCUIT_FAILURE_THRESHOLD:
            return False

        now = time.monotonic()
        holdoff = max(0.0, float(holdoff))
        state.open = True
        state.opened_until = now + holdoff
        state.last_holdoff = holdoff
        state.failure_count = 0
        state.reset_pending = False
        self.last_opened_at = now
        self.last_holdoff = holdoff
        self.last_opened_endpoint = endpoint
        return True


_GTOP_CIRCUIT_FAILURE_THRESHOLD = 5
_GTOP_CIRCUIT_HOLDOFF_SECONDS = 600.0
_GTOP_CIRCUIT_BREAKER = _CircuitBreaker()
_GTOP_CIRCUIT_SKIP_LOG: set[tuple[str, str]] = set()


def _reset_gtop_circuit_state() -> None:
    """Utility for tests to reset the circuit breaker and skip log."""

    _GTOP_CIRCUIT_BREAKER.reset()
    _GTOP_CIRCUIT_SKIP_LOG.clear()


def _collect_name_fields(name_obj: dict[str, Any]) -> Iterable[str]:
    """Yield all full and short names from a UniProt name object."""
    if not isinstance(name_obj, dict):
        return []
    names: list[str] = []
    full = name_obj.get("fullName")
    if isinstance(full, dict):
        value = full.get("value")
        if value:
            names.append(value)
    short = name_obj.get("shortName") or name_obj.get("shortNames")
    if isinstance(short, dict):
        value = short.get("value")
        if value:
            names.append(value)
    elif isinstance(short, list):
        for item in short:
            if isinstance(item, dict):
                value = item.get("value")
                if value:
                    names.append(value)
    return names


def _extract_protein_names(desc: dict[str, Any]) -> set[str]:
    names: set[str] = set()
    if not isinstance(desc, dict):
        return names
    rec = desc.get("recommendedName")
    if isinstance(rec, dict):
        names.update(_collect_name_fields(rec))
    for key in ("alternativeNames", "submissionNames", "submittedName"):
        items = desc.get(key) or []
        for item in items:
            names.update(_collect_name_fields(item))
    return names


def _extract_gene_names(entry: dict[str, Any]) -> set[str]:
    names: set[str] = set()
    for gene in entry.get("genes", []):
        if not isinstance(gene, dict):
            continue
        main = gene.get("geneName")
        if isinstance(main, dict):
            value = main.get("value")
            if value:
                names.add(value)
        for syn in gene.get("synonyms", []):
            if isinstance(syn, dict):
                value = syn.get("value")
                if value:
                    names.add(value)
    return names


def extract_names(data: Any) -> set[str]:
    """Return all protein and gene names found in ``data``.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        A set of name strings aggregated from protein and gene sections.

    """
    names: set[str] = set()
    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        names.update(_extract_protein_names(entry.get("proteinDescription", {})))
        names.update(_extract_gene_names(entry))
    return names


def extract_organism(data: Any) -> dict[str, str]:
    """Return organism taxonomy information for the entry in ``data``.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        A dictionary with keys ``genus``, ``superkingdom``, ``phylum``,
        ``lineage_class`` and ``taxon_id``. Empty strings are returned when a
        field is missing.

    """
    result = {
        "genus": "",
        "superkingdom": "",
        "phylum": "",
        "lineage_class": "",
        "taxon_id": "",
    }
    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        org = entry.get("organism", {})
        if not isinstance(org, dict):
            continue
        taxon_id = org.get("taxonId")
        if taxon_id is not None:
            result["taxon_id"] = str(taxon_id)
        lineage = org.get("lineage") or []
        if isinstance(lineage, list) and lineage:
            result["superkingdom"] = lineage[0]
            phylum_idx = 1
            if len(lineage) >= 2:
                candidate = lineage[1]
                if (
                    isinstance(candidate, str)
                    and candidate.endswith("zoa")
                    and len(lineage) >= 3
                ):
                    result["phylum"] = lineage[2]
                    phylum_idx = 2
                else:
                    result["phylum"] = candidate if isinstance(candidate, str) else ""
            class_idx = phylum_idx + 1
            if class_idx < len(lineage):
                candidate_class = lineage[class_idx]
                if isinstance(candidate_class, str):
                    result["lineage_class"] = candidate_class
            result["genus"] = lineage[-1]
        sci_name = org.get("scientificName")
        if sci_name and not result["genus"]:
            result["genus"] = sci_name.split()[0]
        break
    return result


def _first_entry(data: Any) -> dict[str, Any] | None:
    """Return the first UniProt entry from *data* when available."""

    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if isinstance(entry, dict):
            return entry
    return None


def extract_uniprotkb_id(
    data: Any, *, allow_primary: bool = False, default: str | None = None
) -> str | None:
    """Return the ``uniProtkbId`` (optionally with fallbacks) for ``data``.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.
        allow_primary: When ``True`` and ``uniProtkbId`` is missing, fall back to
            ``primaryAccession``.
        default: Value returned when neither ``uniProtkbId`` nor
            ``primaryAccession`` (when ``allow_primary`` is enabled) is present.

    Returns:
        The ``uniProtkbId`` string when present. If ``allow_primary`` is set,
        the ``primaryAccession`` value is returned when ``uniProtkbId`` is
        absent. ``default`` is used as a final fallback.

    """

    entry = _first_entry(data)
    if not entry:
        return default
    value = entry.get("uniProtkbId")
    if isinstance(value, str) and value:
        return value
    if allow_primary:
        primary = entry.get("primaryAccession")
        if isinstance(primary, str) and primary:
            return primary
    return default


def extract_secondary_accessions(data: Any) -> list[str]:
    """Return secondary accession IDs from ``data``.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        A sorted list of secondary accession identifiers. An empty list is
        returned when no secondary accessions are present.

    """
    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if isinstance(entry, dict):
            secs = entry.get("secondaryAccessions") or []
            if isinstance(secs, list):
                return sorted([s for s in secs if isinstance(s, str)])
            break
    return []


def extract_recommended_name(data: Any) -> str | None:
    """Return the recommended protein name from ``data``.

    The recommended name is taken from
    ``proteinDescription.recommendedName.fullName.value`` when available. If
    that field is missing, the first short name is used instead.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        The recommended name string, or ``None`` when unavailable.

    """
    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        desc = entry.get("proteinDescription", {})
        if not isinstance(desc, dict):
            break
        rec = desc.get("recommendedName")
        if not isinstance(rec, dict):
            break
        full = rec.get("fullName")
        if isinstance(full, dict):
            value = full.get("value")
            if isinstance(value, str):
                return value
        shorts = rec.get("shortNames") or rec.get("shortName")
        if isinstance(shorts, list):
            for item in shorts:
                if isinstance(item, dict):
                    value = item.get("value")
                    if isinstance(value, str):
                        return value
        elif isinstance(shorts, dict):
            value = shorts.get("value")
            if isinstance(value, str):
                return value
        break
    return None


def extract_gene_name(data: Any) -> str | None:
    """Return the primary gene name from ``data``.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        The primary gene name string, or ``None`` when unavailable.

    """
    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        genes = entry.get("genes", [])
        if not isinstance(genes, list):
            break
        for gene in genes:
            if not isinstance(gene, dict):
                continue
            gname = gene.get("geneName")
            if isinstance(gname, dict):
                value = gname.get("value")
                if isinstance(value, str):
                    return value
        break
    return None


def extract_names_for_secondary_accessions(data: Any, *, cfg: UniprotCfg) -> str:
    """Return protein names for secondary accessions listed in ``data``.

    The function looks up each secondary accession via :func:`fetch_uniprot`
    and aggregates all protein names from common description fields. Names are
    deduplicated and returned as a single pipe-separated string. When no names
    are found or the entry cannot be retrieved, an empty string is returned.

    Parameters
    ----------
    data:
        A UniProt JSON structure, list of entries, or search results containing
        UniProt entries.

    Returns
    -------
    str
        Pipe-separated protein names for all secondary accessions.

    """
    names: set[str] = set()
    for acc in extract_secondary_accessions(data):
        try:
            entry = fetch_uniprot(acc, cfg=cfg)
        except UniProtFetchError as exc:  # pragma: no cover - network errors
            logger.warning(
                "secondary_accession_fetch_failed", accession=acc, error=str(exc)
            )
            continue
        desc = entry.get("proteinDescription")
        if isinstance(desc, dict):
            names.update(_extract_protein_names(desc))
    return "|".join(sorted(names))


def _collect_ec_numbers(name_obj: dict[str, Any]) -> Iterable[str]:
    """Yield EC numbers from a UniProt name object."""
    if not isinstance(name_obj, dict):
        return []
    numbers: list[str] = []
    ec = name_obj.get("ecNumbers") or name_obj.get("ecNumber")
    if isinstance(ec, list):
        for item in ec:
            if isinstance(item, dict):
                value = item.get("value")
                if value:
                    numbers.append(value)
            elif isinstance(item, str):
                numbers.append(item)
    elif isinstance(ec, dict):
        value = ec.get("value")
        if value:
            numbers.append(value)
    elif isinstance(ec, str):
        numbers.append(ec)
    return numbers


def extract_keywords(data: Any) -> dict[str, Any]:
    """Return keyword and feature information found in ``data``.

    The function gathers functional keywords, EC numbers, subcellular
    locations, topology hints, and whether transmembrane or intramembrane
    regions are annotated for the entry.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        A dictionary with keys ``molecular_function``, ``cellular_component``,
        ``ec_numbers``, ``subcellular_location``, ``topology``,
        ``transmembrane``, and ``intramembrane``. Keyword-related values are
        returned as sets, while ``transmembrane`` and ``intramembrane`` are
        booleans.

    """
    result: dict[str, Any] = {
        "molecular_function": set(),
        "cellular_component": set(),
        "ec_numbers": set(),
        "subcellular_location": set(),
        "topology": set(),
        "transmembrane": False,
        "intramembrane": False,
    }
    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if not isinstance(entry, dict):
            continue

        # Keyword extraction
        for kw in entry.get("keywords", []):
            if not isinstance(kw, dict):
                continue
            category = kw.get("category")
            if isinstance(category, dict):
                category = category.get("value")
            name = kw.get("name")
            if isinstance(name, dict):
                name = name.get("value")
            if not isinstance(name, str):
                continue
            if category == "Molecular function":
                result["molecular_function"].add(name)
            elif category == "Cellular component":
                result["cellular_component"].add(name)

        # EC number extraction from protein descriptions
        desc = entry.get("proteinDescription", {})
        if isinstance(desc, dict):
            rec = desc.get("recommendedName")
            if isinstance(rec, dict):
                result["ec_numbers"].update(_collect_ec_numbers(rec))
            for key in ("alternativeNames", "submissionNames"):
                items = desc.get(key) or []
                for item in items:
                    result["ec_numbers"].update(_collect_ec_numbers(item))

        # Subcellular location and topology
        comments = entry.get("comments", [])
        if isinstance(comments, list):
            for comment in comments:
                if not isinstance(comment, dict):
                    continue
                if comment.get("commentType") != "SUBCELLULAR LOCATION":
                    continue
                sublocs = comment.get("subcellularLocations") or []
                for loc in sublocs:
                    if not isinstance(loc, dict):
                        continue
                    sub = loc.get("location")
                    if isinstance(sub, dict):
                        value = sub.get("value")
                        if isinstance(value, str):
                            result["subcellular_location"].add(value)
                    topo = loc.get("topology")
                    if isinstance(topo, dict):
                        value = topo.get("value")
                        if isinstance(value, str):
                            result["topology"].add(value)

        # Feature flags for membranes
        features = entry.get("features", [])
        if isinstance(features, list):
            for feat in features:
                if not isinstance(feat, dict):
                    continue
                ftype = feat.get("type")
                if ftype == "TRANSMEMBRANE":
                    result["transmembrane"] = True
                elif ftype == "INTRAMEMBRANE":
                    result["intramembrane"] = True
    return result


def extract_ptm(data: Any) -> dict[str, bool]:
    """Return post-translational modification flags found in ``data``.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        A dictionary mapping PTM feature names to booleans. The returned keys
        include ``glycosylation``, ``lipidation``, ``disulfide_bond``,
        ``modified_residue``, ``phosphorylation``, ``acetylation``,
        ``ubiquitination``, ``signal_peptide``, ``propeptide``, and
        ``transmembrane``.

    """
    feature_map = {
        "glycosylation": "GLYCOSYLATION",
        "lipidation": "LIPIDATION",
        "disulfide_bond": "DISULFIDE BOND",
        "modified_residue": "MODIFIED RESIDUE",
        "phosphorylation": "PHOSPHORYLATION",
        "acetylation": "ACETYLATION",
        "ubiquitination": "UBIQUITINATION",
        "signal_peptide": "SIGNAL PEPTIDE",
        "propeptide": "PROPEPTIDE",
        "transmembrane": "TRANSMEMBRANE",
    }
    result: dict[str, bool] = {key: False for key in feature_map}
    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        features = entry.get("features", [])
        if not isinstance(features, list):
            continue
        for feat in features:
            if not isinstance(feat, dict):
                continue
            ftype = feat.get("type")
            if not isinstance(ftype, str):
                continue
            uftype = ftype.upper()
            for key, expected in feature_map.items():
                if uftype == expected.upper():
                    result[key] = True
    return result


def extract_isoform(data: Any) -> dict[str, str]:
    """Return isoform information found in ``data``.

    The function inspects ``ALTERNATIVE PRODUCTS`` comments and gathers the
    names, IDs, and synonyms for each isoform. Multiple IDs or synonyms within
    an isoform are joined by ``":"`` while separate isoforms are joined by
    ``"|"``. When no isoform data is available, the strings ``"None"`` are
    returned for all fields.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        A dictionary with keys ``isoform_names``, ``isoform_ids``, and
        ``isoform_synonyms`` mapping to pipe separated strings.

    """
    names: list[str] = []
    ids: list[str] = []
    syns: list[str] = []
    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        comments = entry.get("comments", [])
        if not isinstance(comments, list):
            continue
        for comment in comments:
            if (
                not isinstance(comment, dict)
                or comment.get("commentType") != "ALTERNATIVE PRODUCTS"
            ):
                continue
            isoforms = comment.get("isoforms", [])
            if not isinstance(isoforms, list):
                continue
            for iso in isoforms:
                if not isinstance(iso, dict):
                    continue
                name = None
                name_obj = iso.get("name")
                if isinstance(name_obj, dict):
                    name = name_obj.get("value")
                if isinstance(name, str):
                    names.append(name)
                # IDs
                iso_ids: list[str] = []
                for iid in iso.get("isoformIds", []) or []:
                    if isinstance(iid, str):
                        iso_ids.append(iid)
                ids.append(":".join(iso_ids) if iso_ids else "N/A")
                # Synonyms
                syn_list: list[str] = []
                for syn in iso.get("synonyms", []) or []:
                    if isinstance(syn, dict):
                        value = syn.get("value")
                        if isinstance(value, str):
                            syn_list.append(value)
                syns.append(":".join(syn_list) if syn_list else "N/A")
    result = {
        "isoform_names": "|".join(names) if names else "None",
        "isoform_ids": "|".join(ids) if names else "None",
        "isoform_synonyms": "|".join(syns) if names else "None",
    }
    return result


def extract_crossrefs(data: Any) -> dict[str, str]:
    """Return cross-reference identifiers for selected databases.

    The UniProt record contains a list of cross references for many external
    databases. This helper searches that list and aggregates the identifiers for
    a predefined subset of databases. When multiple identifiers are present for
    the same database, they are returned as a pipe-separated string. Missing
    databases yield empty strings.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        A dictionary mapping database names to pipe-separated identifier
        strings. The returned keys include ``GuidetoPHARMACOLOGY``, ``family``,
        ``SUPFAM``, ``PROSITE``, ``InterPro``, ``Pfam``, ``PRINTS``, ``TCDB``,
        ``xref_pdb``, ``xref_alphafold``, and ``xref_ensembl``.

    """
    db_map: dict[str, tuple[str, ...]] = {
        "GuidetoPHARMACOLOGY": ("GuidetoPHARMACOLOGY",),
        "family": ("family",),
        "SUPFAM": ("SUPFAM",),
        "PROSITE": ("PROSITE",),
        "InterPro": ("InterPro",),
        "Pfam": ("Pfam",),
        "PRINTS": ("PRINTS",),
        "TCDB": ("TCDB",),
        "xref_pdb": ("PDB",),
        "xref_alphafold": ("AlphaFoldDB",),
        "xref_ensembl": ("Ensembl",),
    }
    result: dict[str, list[str]] = {column: [] for column in db_map}
    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        xrefs = (
            entry.get("uniProtKBCrossReferences")
            or entry.get("uniProtCrossReferences")
            or entry.get("dbReferences")
            or []
        )
        if not isinstance(xrefs, list):
            continue
        for ref in xrefs:
            if not isinstance(ref, dict):
                continue
            db = ref.get("database")
            if not isinstance(db, str):
                continue
            ref_id = ref.get("id")
            if not isinstance(ref_id, str):
                continue
            for column, names in db_map.items():
                if db in names:
                    result[column].append(ref_id)
                    break
    return {column: "|".join(ids) for column, ids in result.items()}


def extract_activity(data: Any) -> dict[str, str]:
    """Return catalytic reaction names and EC numbers found in ``data``.

    The UniProt record may list one or more "CATALYTIC ACTIVITY" comments,
    each describing a reaction and an associated EC number. This helper
    aggregates those reactions and numbers as pipe-separated strings.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        A dictionary with keys ``reactions`` and ``reaction_ec_numbers``.
        Missing information yields empty strings.

    """
    reactions: list[str] = []
    numbers: list[str] = []
    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        comments = entry.get("comments", [])
        if not isinstance(comments, list):
            continue
        for comment in comments:
            if not isinstance(comment, dict):
                continue
            if comment.get("commentType") != "CATALYTIC ACTIVITY":
                continue
            reaction = comment.get("reaction")
            if not isinstance(reaction, dict):
                continue
            name = reaction.get("name")
            if isinstance(name, dict):
                name = name.get("value")
            if isinstance(name, str):
                reactions.append(name)
            numbers.extend(list(_collect_ec_numbers(reaction)))
    return {
        "reactions": "|".join(reactions),
        "reaction_ec_numbers": "|".join(numbers),
    }


def _fetch_gtop_endpoint(
    gtop_id: str,
    endpoint: str,
    *,
    cfg: IupharCfg,
) -> Any:
    """Return JSON payload for ``endpoint`` of a Guide-to-Pharmacology target."""

    if not getattr(cfg, "enable", True):
        logger.debug(
            "gtop_fetch_disabled",
            gtop_id=gtop_id,
            endpoint=endpoint,
        )
        return []

    cache_key = (gtop_id, endpoint)
    if cache_key in _GTOP_JSON_FAILURE_CACHE:
        if cache_key not in _GTOP_SKIPPED_FAILURE_LOG:
            logger.debug(
                "gtop_fetch_skipped_after_failure",
                gtop_id=gtop_id,
                endpoint=endpoint,
            )
            _GTOP_SKIPPED_FAILURE_LOG.add(cache_key)
        return None

    decision = _GTOP_CIRCUIT_BREAKER.before_call(endpoint)
    if decision.reset_holdoff is not None:
        if _GTOP_CIRCUIT_SKIP_LOG:
            _GTOP_CIRCUIT_SKIP_LOG.clear()
        logger.info(
            "gtop_circuit_reset",
            downtime=decision.reset_holdoff,
        )

    if not decision.allow_call:
        if cache_key not in _GTOP_CIRCUIT_SKIP_LOG:
            logger.warning(
                "gtop_circuit_open_skip",
                gtop_id=gtop_id,
                endpoint=endpoint,
                retry_after=round(decision.remaining, 3),
            )
            _GTOP_CIRCUIT_SKIP_LOG.add(cache_key)
        return None

    limiter = get_limiter("iuphar", cfg.rps, cfg.burst)
    base = cfg.base.rstrip("/")
    path = f"/{endpoint.lstrip('/')}" if endpoint else ""
    url = f"{base}/targets/{gtop_id}{path}"
    timeout = (cfg.timeout_connect, cfg.timeout_read)
    limiter.acquire()
    try:
        session = get_uniprot_session()
        with session.get(url, timeout=timeout) as response:
            status_code = getattr(response, "status_code", None)
            if status_code == 404:
                logger.info(
                    "gtop_endpoint_missing",
                    gtop_id=gtop_id,
                    endpoint=endpoint,
                    status_code=status_code,
                )
                _GTOP_JSON_FAILURE_CACHE.add(cache_key)
                _GTOP_CIRCUIT_BREAKER.record_success(endpoint)
                return None

            response.raise_for_status()
            raw_content_type = response.headers.get("Content-Type")
            content_type = raw_content_type if isinstance(raw_content_type, str) else ""
            body = response.content
            if isinstance(body, bytes | bytearray):
                body_is_empty = not body.strip()
            elif isinstance(body, str):
                body_is_empty = not body.strip()
            else:
                body_is_empty = not body
            content_length = response.headers.get("Content-Length")
            if content_length == "0" or body_is_empty:
                logger.info(
                    "gtop_empty_response",
                    gtop_id=gtop_id,
                    endpoint=endpoint,
                    content_type=raw_content_type,
                )
                _GTOP_CIRCUIT_BREAKER.record_success(endpoint)
                return []
            try:
                payload = response.json()
            except (json.JSONDecodeError, ValueError) as exc:
                logger.warning(
                    "gtop_json_decode_failed",
                    gtop_id=gtop_id,
                    endpoint=endpoint,
                    error=str(exc),
                    content_type=raw_content_type,
                )
                _GTOP_JSON_FAILURE_CACHE.add(cache_key)
                opened = _GTOP_CIRCUIT_BREAKER.record_failure(
                    endpoint, _GTOP_CIRCUIT_HOLDOFF_SECONDS
                )
                if opened:
                    logger.warning(
                        "gtop_circuit_opened",
                        gtop_id=gtop_id,
                        endpoint=endpoint,
                        retry_after=_GTOP_CIRCUIT_BREAKER.last_holdoff,
                        failure_threshold=_GTOP_CIRCUIT_FAILURE_THRESHOLD,
                    )
                return None
            if "json" not in content_type.lower():
                if cache_key not in _GTOP_NON_JSON_CONTENT_TYPE_CACHE:
                    logger.warning(
                        "gtop_non_json_response",
                        gtop_id=gtop_id,
                        endpoint=endpoint,
                        content_type=raw_content_type,
                    )
                    _GTOP_NON_JSON_CONTENT_TYPE_CACHE.add(cache_key)
            _GTOP_CIRCUIT_BREAKER.record_success(endpoint)
            return payload
    except requests.RequestException as exc:  # pragma: no cover - network failures
        response = getattr(exc, "response", None)
        status_code = getattr(response, "status_code", None)
        if status_code == 404:
            logger.info(
                "gtop_endpoint_missing",
                gtop_id=gtop_id,
                endpoint=endpoint,
                status_code=status_code,
            )
            _GTOP_JSON_FAILURE_CACHE.add(cache_key)
            _GTOP_CIRCUIT_BREAKER.record_success(endpoint)
            return None
        logger.warning(
            "gtop_request_failed",
            gtop_id=gtop_id,
            endpoint=endpoint,
            error=str(exc),
            status_code=status_code,
        )
        if status_code is None or status_code >= 500:
            opened = _GTOP_CIRCUIT_BREAKER.record_failure(
                endpoint, _GTOP_CIRCUIT_HOLDOFF_SECONDS
            )
            if opened:
                logger.warning(
                    "gtop_circuit_opened",
                    gtop_id=gtop_id,
                    endpoint=endpoint,
                    retry_after=_GTOP_CIRCUIT_BREAKER.last_holdoff,
                    failure_threshold=_GTOP_CIRCUIT_FAILURE_THRESHOLD,
                )
        else:
            _GTOP_CIRCUIT_BREAKER.record_success(endpoint)
        _GTOP_JSON_FAILURE_CACHE.add(cache_key)
    return None


def _summarise_gtop_function(entries: Any) -> str:
    """Return a concise textual summary from a function payload."""

    if not isinstance(entries, list):
        return ""
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        description = entry.get("description")
        property_text = entry.get("property")
        if isinstance(property_text, str) and property_text.strip():
            if isinstance(description, str) and description.strip():
                return f"{description.strip()}: {property_text.strip()}"
            return property_text.strip()
        if isinstance(description, str) and description.strip():
            return description.strip()
        tissue = entry.get("tissue")
        if isinstance(tissue, str) and tissue.strip():
            return tissue.strip()
    return ""


def _update_gtop_metadata(
    result: dict[str, Any],
    *,
    cfg: IupharCfg | None = None,
) -> None:
    """Populate Guide-to-Pharmacology statistics in ``result`` when available."""

    gtop_value = result.get("GuidetoPHARMACOLOGY")
    if not isinstance(gtop_value, str) or not gtop_value:
        return
    gtop_id = gtop_value.split("|", 1)[0].strip()
    if not gtop_id:
        return
    config = cfg or IupharCfg()
    if not getattr(config, "enable", True):
        logger.debug(
            "gtop_enrichment_disabled",
            gtop_id=gtop_id,
        )
        return

    natural = _fetch_gtop_endpoint(gtop_id, "naturalLigands", cfg=config)
    if isinstance(natural, list):
        result["gtop_natural_ligands_n"] = str(len(natural))

    interactions = _fetch_gtop_endpoint(gtop_id, "interactions", cfg=config)
    if isinstance(interactions, list):
        result["gtop_interactions_n"] = str(len(interactions))

    function_entries = _fetch_gtop_endpoint(gtop_id, "function", cfg=config)
    summary = _summarise_gtop_function(function_entries)
    if summary:
        result["gtop_function_text_short"] = summary


def iter_ids(csv_path: str, sep: str = ",", encoding: str = "utf-8") -> Iterable[str]:
    """Yield UniProt IDs from a CSV file with a ``uniprot_id`` column.

    Parameters
    ----------
    csv_path:
        Path to a CSV file containing a ``uniprot_id`` column.
    sep:
        Field delimiter used in ``csv_path``. Defaults to a comma.
    encoding:
        Text encoding of ``csv_path``. Defaults to UTF-8.

    Yields
    ------
    str
        Each UniProt accession ID.

    """
    try:
        with open(csv_path, newline="", encoding=encoding) as handle:
            reader = csv.DictReader(handle, delimiter=sep)
            if reader.fieldnames is None or "uniprot_id" not in reader.fieldnames:
                raise ValueError("Input CSV must have a uniprot_id column")
            for row in reader:
                uid = row.get("uniprot_id")
                if uid:
                    yield uid.strip()
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"input file not found: {csv_path}") from exc
    except csv.Error as exc:
        raise ValueError(f"malformed CSV in file: {csv_path}: {exc}") from exc


def _extract_audit_last_update(audit: dict[str, Any]) -> str | None:
    """Return the latest available audit date for a UniProt entry.

    UniProt's JSON payloads may expose the audit dates using several field
    names depending on the release. ``lastUpdateDate`` is preferred, but older
    or transitional payloads can include ``lastAnnotationUpdateDate`` or
    ``lastSequenceUpdateDate``. Each of those fields can be either a plain
    string or an object containing a ``value`` key. This helper returns the
    first populated date found in that priority order.

    Parameters
    ----------
    audit:
        ``entryAudit`` section from a UniProt JSON entry.

    Returns
    -------
    str | None
        The best available audit date, or ``None`` if all candidates are
        missing.
    """

    def _coerce_date(value: Any) -> str | None:
        if isinstance(value, str):
            return value
        if isinstance(value, dict):
            candidate = value.get("value")
            if isinstance(candidate, str):
                return candidate
        return None

    for field_name in (
        "lastUpdateDate",
        "lastAnnotationUpdateDate",
        "lastSequenceUpdateDate",
    ):
        date_value = _coerce_date(audit.get(field_name))
        if date_value:
            return date_value
    return None


def _load_uniprot_entry(
    uid: str, data_dir: Path, cfg: UniprotCfg
) -> tuple[Any | None, str | None]:
    """Return UniProt JSON data for ``uid`` with isoform-aware fallback."""

    fallback_id = _canonical_accession(uid)
    candidates: list[str] = [uid]
    if fallback_id and fallback_id not in candidates:
        candidates.append(fallback_id)

    for candidate in candidates:
        json_path = data_dir / f"{candidate}.json"
        try:
            with open(json_path, encoding="utf-8") as handle:
                data = json.load(handle)
        except FileNotFoundError:
            continue
        except json.JSONDecodeError:
            logger.warning("uniprot_json_malformed", uid=candidate)
            continue
        else:
            if candidate != uid:
                logger.info(
                    "uniprot_isoform_fallback_cache",
                    isoform=uid,
                    canonical=candidate,
                )
            return data, candidate

    logger.info("uniprot_json_download", uid=uid)
    for index, candidate in enumerate(candidates):
        try:
            data = fetch_uniprot(candidate, cfg=cfg)
        except UniProtFetchError as exc:
            if candidate == uid and fallback_id and index == 0:
                logger.info(
                    "uniprot_isoform_fetch_retry",
                    isoform=uid,
                    canonical=fallback_id,
                    error=str(exc),
                )
                continue
            logger.warning("uniprot_json_fetch_failed", uid=candidate, error=str(exc))
            return None, None

        data_dir.mkdir(parents=True, exist_ok=True)
        persist_ids = [candidate]
        if candidate != uid:
            persist_ids.append(uid)
            logger.info(
                "uniprot_isoform_fallback_fetch",
                isoform=uid,
                canonical=candidate,
            )
        for persist_id in persist_ids:
            json_path = data_dir / f"{persist_id}.json"
            try:
                with open(json_path, "w", encoding="utf-8") as handle:
                    json.dump(data, handle)
            except OSError as exc:  # pragma: no cover - disk I/O failure
                logger.warning(
                    "uniprot_json_write_failed", uid=persist_id, error=str(exc)
                )
        return data, candidate

    return None, None


def collect_info(
    uid: str,
    data_dir: Path | str | None = None,
    *,
    cfg: UniprotCfg,
    gtop_cfg: IupharCfg | None = None,
) -> dict[str, Any]:
    """Return names, organism, keyword, PTM, isoform, cross-ref, and activity data for ``uid``.

    Parameters
    ----------
    uid:
        UniProt accession identifier.
    data_dir:
        Directory containing ``<uid>.json`` files with UniProt data. If not
        provided, :data:`_DEFAULT_UNIPROT_DATA_DIR` is used.
    cfg:
        UniProt configuration used for downloading missing records.
    gtop_cfg:
        Guide-to-Pharmacology configuration for enriching cross references.

    Returns
    -------
    dict
        A dictionary with keys ``uniprot_id``, ``names``, organism taxonomy
        fields, keyword categories, EC numbers, subcellular location data,
        membrane features, post-translational modification flags, isoform
        metadata, selected database cross references, and Guide-to-
        Pharmacology statistics. Missing or invalid files leave fields empty.

    """
    if data_dir is None:
        data_dir = _DEFAULT_UNIPROT_DATA_DIR
    data_dir = Path(data_dir)

    result = {
        "uniprot_id": uid,
        "names": "",
        "genus": "",
        "superkingdom": "",
        "phylum": "",
        "lineage_class": "",
        "taxon_id": "",
        "sequence_length": "",
        "molecular_function": "",
        "cellular_component": "",
        "ec_numbers": "",
        "subcellular_location": "",
        "topology": "",
        "transmembrane": False,
        "intramembrane": False,
        "glycosylation": False,
        "lipidation": False,
        "disulfide_bond": False,
        "modified_residue": False,
        "phosphorylation": False,
        "acetylation": False,
        "ubiquitination": False,
        "signal_peptide": False,
        "propeptide": False,
        "isoform_names": "None",
        "isoform_ids": "None",
        "isoform_synonyms": "None",
        "GuidetoPHARMACOLOGY": "",
        "family": "",
        "SUPFAM": "",
        "PROSITE": "",
        "InterPro": "",
        "Pfam": "",
        "PRINTS": "",
        "TCDB": "",
        "reactions": "",
        "reaction_ec_numbers": "",
        "secondaryAccessionNames": "",
        "gtop_natural_ligands_n": "",
        "gtop_interactions_n": "",
        "gtop_function_text_short": "",
        "xref_pdb": "",
        "xref_alphafold": "",
        "xref_ensembl": "",
        "uniprot_last_update": "",
        "uniprot_version": "",
        "pipeline_version": "",
        "timestamp_utc": "",
        "uniProtkbIdFallback": "",
    }
    data, _source_id = _load_uniprot_entry(uid, data_dir, cfg)
    if data is None:
        return result

    names = extract_names(data)
    org = extract_organism(data)
    keywords = extract_keywords(data)
    ptm = extract_ptm(data)
    iso = extract_isoform(data)
    cross = extract_crossrefs(data)
    activity = extract_activity(data)
    entry = data
    if isinstance(entry, dict) and "results" in entry:
        results = entry.get("results")
        entry = results[0] if isinstance(results, list) and results else entry
    elif isinstance(entry, list) and entry:
        entry = entry[0]
    if isinstance(entry, dict):
        sequence = entry.get("sequence")
        if isinstance(sequence, dict):
            length = sequence.get("length")
            if length is not None:
                result["sequence_length"] = str(length)
        audit = entry.get("entryAudit")
        if isinstance(audit, dict):
            last_update = _extract_audit_last_update(audit)
            if last_update:
                result["uniprot_last_update"] = last_update
            version = audit.get("entryVersion")
            if version is not None:
                result["uniprot_version"] = str(version)
    result["names"] = "|".join(sorted(names))
    result.update(org)
    result["molecular_function"] = "|".join(sorted(keywords["molecular_function"]))
    result["cellular_component"] = "|".join(sorted(keywords["cellular_component"]))
    result["ec_numbers"] = "|".join(sorted(keywords["ec_numbers"]))
    result["subcellular_location"] = (
        "|".join(sorted(keywords["subcellular_location"])) or "N/A"
    )
    result["topology"] = "|".join(sorted(keywords["topology"])) or "N/A"
    result["transmembrane"] = ptm["transmembrane"]
    result["intramembrane"] = keywords["intramembrane"]
    for key in (
        "glycosylation",
        "lipidation",
        "disulfide_bond",
        "modified_residue",
        "phosphorylation",
        "acetylation",
        "ubiquitination",
        "signal_peptide",
        "propeptide",
    ):
        result[key] = ptm[key]
    result.update(iso)
    result.update(cross)
    _update_gtop_metadata(result, cfg=gtop_cfg)
    result.update(activity)
    uni_protkb_id = extract_uniprotkb_id(data)
    fallback_uniprotkb_id = uni_protkb_id or extract_uniprotkb_id(
        data, allow_primary=True
    )
    result["uniProtkbId"] = uni_protkb_id
    result["uniProtkbIdFallback"] = fallback_uniprotkb_id or uid
    result["secondaryAccessions"] = extract_secondary_accessions(data)
    result["recommendedName"] = extract_recommended_name(data)
    result["geneName"] = extract_gene_name(data)
    result["secondaryAccessionNames"] = extract_names_for_secondary_accessions(
        data, cfg=cfg
    )
    return result


def process(
    input_csv: str,
    output_csv: str,
    data_dir: Path | str | None = None,
    *,
    cfg: UniprotCfg,
    gtop_cfg: IupharCfg | None = None,
    sep: str = ",",
    encoding: str = "utf-8",
) -> None:
    """Read IDs from ``input_csv`` and write extracted data to ``output_csv``.

    The output includes names, taxonomy, keyword categories, EC numbers,
    subcellular locations, membrane and PTM flags, isoform metadata,
    catalytic reactions with EC numbers, and selected database cross references
    for each accession.

    Parameters
    ----------
    input_csv:
        Path to the CSV file listing UniProt IDs.
    output_csv:
        Destination path for the output CSV file.
    data_dir:
        Directory where JSON files for each ID are stored. Defaults to
        :data:`_DEFAULT_UNIPROT_DATA_DIR`.
    cfg:
        UniProt configuration used for network requests when local files are
        missing.
    gtop_cfg:
        Guide-to-Pharmacology configuration applied when enriching
        cross-reference data.
    sep:
        Field delimiter used for both input and output CSV files. Defaults to a comma.
    encoding:
        File encoding for both input and output CSV files. Defaults to UTF-8.

    Returns
    -------
    None
        The processed information is written to ``output_csv``.

    """
    if data_dir is None:
        data_dir = _DEFAULT_UNIPROT_DATA_DIR
    data_dir = Path(data_dir)

    fieldnames = UNIPROT_OUTPUT_COLUMNS
    expected_columns = set(fieldnames)

    try:
        with open(output_csv, "w", newline="", encoding=encoding) as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=sep)
            writer.writeheader()
            for uid in iter_ids(input_csv, sep=sep, encoding=encoding):
                info = collect_info(uid, data_dir, cfg=cfg, gtop_cfg=gtop_cfg)
                unexpected = sorted(set(info) - expected_columns)
                if unexpected:
                    logger.debug("uniprot_extra_fields", uid=uid, columns=unexpected)
                row = {
                    column: ("" if info.get(column) is None else info.get(column))
                    for column in fieldnames
                }
                secondary = row.get("secondaryAccessions")
                if isinstance(secondary, list):
                    row["secondaryAccessions"] = "|".join(secondary)
                writer.writerow(row)
    except OSError as exc:
        raise OSError(f"failed to write output CSV: {output_csv}: {exc}") from exc
