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
    Given a UniProt accession and directory containing ``<uid>.json``
    files, return a dictionary with the accession, all names, and
    organism taxonomy data.

``process(input_csv, output_csv, data_dir="uniprot")``
    Batch-process a CSV of UniProt IDs and write an output CSV with
    names and organism information for each ID.
"""

from __future__ import annotations

import csv
import json
import logging
import os
from typing import Any, Dict, Iterable, List, Set

import requests
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)

# Shared HTTP session with retry/backoff to make network calls more robust.
_retry = Retry(
    total=3,
    backoff_factor=1.0,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET"],
)
_session: Session = requests.Session()
_session.mount("http://", HTTPAdapter(max_retries=_retry))
_session.mount("https://", HTTPAdapter(max_retries=_retry))

API_URL = "https://rest.uniprot.org/uniprotkb/{id}.json"

__all__ = [
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
]


def fetch_uniprot(uniprot_id: str) -> Dict[str, Any]:
    """Fetch a UniProt JSON record from the public REST API.

    Parameters
    ----------
    uniprot_id:
        UniProt accession identifier to retrieve.

    Returns
    -------
    dict
        JSON-decoded response.  An empty dictionary is returned when the
        request fails or the payload cannot be decoded.
    """

    url = API_URL.format(id=uniprot_id)
    try:
        resp = _session.get(url, timeout=30)
        resp.raise_for_status()
        try:
            return resp.json()
        except json.JSONDecodeError as exc:  # pragma: no cover - malformed JSON
            logger.warning("Failed to decode JSON for UniProt %s: %s", uniprot_id, exc)
            return {}
    except requests.RequestException as exc:  # pragma: no cover - network
        logger.warning("UniProt request failed for %s: %s", uniprot_id, exc)
        return {}


def _collect_name_fields(name_obj: Dict[str, Any]) -> Iterable[str]:
    """Yield all full and short names from a UniProt name object."""
    if not isinstance(name_obj, dict):
        return []
    names: List[str] = []
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


def _extract_protein_names(desc: Dict[str, Any]) -> Set[str]:
    names: Set[str] = set()
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


def _extract_gene_names(entry: Dict[str, Any]) -> Set[str]:
    names: Set[str] = set()
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


def extract_names(data: Any) -> Set[str]:
    """Return all protein and gene names found in ``data``.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        A set of name strings aggregated from protein and gene sections.
    """
    names: Set[str] = set()
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


def extract_organism(data: Any) -> Dict[str, str]:
    """Return organism taxonomy information for the entry in ``data``.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        A dictionary with keys ``genus``, ``superkingdom``, ``phylum`` and
        ``taxon_id``. Empty strings are returned when a field is missing.
    """
    result = {"genus": "", "superkingdom": "", "phylum": "", "taxon_id": ""}
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
            if len(lineage) >= 2:
                candidate = lineage[1]
                if (
                    isinstance(candidate, str)
                    and candidate.endswith("zoa")
                    and len(lineage) >= 3
                ):
                    result["phylum"] = lineage[2]
                else:
                    result["phylum"] = candidate
            result["genus"] = lineage[-1]
        sci_name = org.get("scientificName")
        if sci_name and not result["genus"]:
            result["genus"] = sci_name.split()[0]
        break
    return result

def extract_uniprotkb_id(data: Any) -> str | None:
    """Return the ``uniProtkbId`` for the first entry in ``data``.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        The ``uniProtkbId`` string when present, otherwise ``None``.
    """

    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if isinstance(entry, dict):
            value = entry.get("uniProtkbId")
            if isinstance(value, str):
                return value
            break
    return None


def extract_secondary_accessions(data: Any) -> List[str]:
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

def extract_names_for_secondary_accessions(data: Any) -> str:
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

    names: Set[str] = set()
    for acc in extract_secondary_accessions(data):
        entry = fetch_uniprot(acc)
        if not isinstance(entry, dict):
            continue
        desc = entry.get("proteinDescription")
        if isinstance(desc, dict):
            names.update(_extract_protein_names(desc))
    return "|".join(sorted(names))

def _collect_ec_numbers(name_obj: Dict[str, Any]) -> Iterable[str]:
    """Yield EC numbers from a UniProt name object."""
    if not isinstance(name_obj, dict):
        return []
    numbers: List[str] = []
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


def extract_keywords(data: Any) -> Dict[str, Any]:
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
    result: Dict[str, Any] = {
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


def extract_ptm(data: Any) -> Dict[str, bool]:
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
    result: Dict[str, bool] = {key: False for key in feature_map}
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


def extract_isoform(data: Any) -> Dict[str, str]:
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
    names: List[str] = []
    ids: List[str] = []
    syns: List[str] = []
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
            if not isinstance(comment, dict) or comment.get("commentType") != "ALTERNATIVE PRODUCTS":
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
                iso_ids: List[str] = []
                for iid in iso.get("isoformIds", []) or []:
                    if isinstance(iid, str):
                        iso_ids.append(iid)
                ids.append(":".join(iso_ids) if iso_ids else "N/A")
                # Synonyms
                syn_list: List[str] = []
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


def extract_crossrefs(data: Any) -> Dict[str, str]:
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
        ``SUPFAM``, ``PROSITE``, ``InterPro``, ``Pfam``, ``PRINTS``, and
        ``TCDB``.
    """
    dbs = [
        "GuidetoPHARMACOLOGY",
        "family",
        "SUPFAM",
        "PROSITE",
        "InterPro",
        "Pfam",
        "PRINTS",
        "TCDB",
    ]
    result: Dict[str, List[str]] = {db: [] for db in dbs}
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
            if db in result:
                ref_id = ref.get("id")
                if isinstance(ref_id, str):
                    result[db].append(ref_id)
    return {db: "|".join(ids) for db, ids in result.items()}


def extract_activity(data: Any) -> Dict[str, str]:
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
    reactions: List[str] = []
    numbers: List[str] = []
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

def collect_info(uid: str, data_dir: str = "uniprot") -> Dict[str, Any]:
    """Return names, organism, keyword, PTM, isoform, cross-ref, and activity data for ``uid``.

    Args:
        uid: UniProt accession identifier.
        data_dir: Directory containing ``<uid>.json`` files with UniProt data.

    Returns:
        A dictionary with keys ``uniprot_id``, ``names``, organism taxonomy
        fields, keyword categories, EC numbers, subcellular location data,
        membrane features, post-translational modification flags, isoform
        metadata, and selected database cross references. Missing or invalid
        files leave fields empty.
    """
    json_path = os.path.join(data_dir, f"{uid}.json")
    result = {
        "uniprot_id": uid,
        "names": "",
        "genus": "",
        "superkingdom": "",
        "phylum": "",
        "taxon_id": "",
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
    }
    try:
        with open(json_path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except FileNotFoundError:
        logger.info("downloading UniProt JSON for %s", uid)
        data = fetch_uniprot(uid)
        if not data:
            logger.warning("failed to retrieve UniProt JSON for %s", uid)
            return result
        os.makedirs(data_dir, exist_ok=True)
        try:
            with open(json_path, "w", encoding="utf-8") as handle:
                json.dump(data, handle)
        except OSError as exc:  # pragma: no cover - disk I/O failure
            logger.warning("unable to write UniProt JSON for %s: %s", uid, exc)
            return result
    except json.JSONDecodeError:
        logger.warning("malformed UniProt JSON for %s", uid)
        return result

    names = extract_names(data)
    org = extract_organism(data)
    keywords = extract_keywords(data)
    ptm = extract_ptm(data)
    iso = extract_isoform(data)
    cross = extract_crossrefs(data)
    activity = extract_activity(data)
    result["names"] = "|".join(sorted(names))
    result.update(org)
    result["molecular_function"] = "|".join(
        sorted(keywords["molecular_function"])
    )
    result["cellular_component"] = "|".join(
        sorted(keywords["cellular_component"])
    )
    result["ec_numbers"] = "|".join(sorted(keywords["ec_numbers"]))
    result["subcellular_location"] = "|".join(
        sorted(keywords["subcellular_location"])
    ) or "N/A"
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
    result.update(activity)
    result["uniProtkbId"] = extract_uniprotkb_id(data)
    result["secondaryAccessions"] = extract_secondary_accessions(data)
    result["recommendedName"] = extract_recommended_name(data)
    result["geneName"] = extract_gene_name(data)
    result["secondaryAccessionNames"] = extract_names_for_secondary_accessions(data)
    return result


def process(
    input_csv: str,
    output_csv: str,
    data_dir: str = "uniprot",
    *,
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
        Directory where JSON files for each ID are stored.
    sep:
        Field delimiter used for both input and output CSV files. Defaults to a comma.
    encoding:
        File encoding for both input and output CSV files. Defaults to UTF-8.

    Returns
    -------
    None
        The processed information is written to ``output_csv``.
    """

    fieldnames = [
        "uniprot_id",
        "names",
        "genus",
        "superkingdom",
        "phylum",
        "taxon_id",
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
        "uniProtkbId",
        "secondaryAccessions",
        "recommendedName",
        "geneName",
        "secondaryAccessionNames",
    ]

    try:
        with open(output_csv, "w", newline="", encoding=encoding) as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=sep)
            writer.writeheader()
            for uid in iter_ids(input_csv, sep=sep, encoding=encoding):
                info = collect_info(uid, data_dir)
                info["secondaryAccessions"] = "|".join(info["secondaryAccessions"])
                writer.writerow(info)
    except OSError as exc:
        raise OSError(f"failed to write output CSV: {output_csv}: {exc}") from exc
