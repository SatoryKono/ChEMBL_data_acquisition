"""
=================

Utility functions to query the ChEMBL API for target and assay
information.  This module rewrites the original Power Query (M)
implementation in Python and exposes a small set of helpers for
retrieving data and augmenting :class:`pandas.DataFrame` objects.

The functions are intentionally resilient: network errors or malformed
responses result in empty structures being returned.  Errors are logged so
that callers may inspect the cause.

The public API is comprised of the following functions:

``get_target``
    Fetch a single target description.
``get_targets``
    Retrieve multiple target descriptions at once.
``get_assay``
    Fetch information about a single assay.
``get_assays``
    Retrieve multiple assay descriptions at once.
``extend_target``
    Augment a DataFrame with columns returned from ``get_target``.

Example
-------
>>> from script.chembl_lib import get_target
>>> get_target("CHEMBL25")
{'pref_name': 'Mu opioid receptor', ...}

The module uses :mod:`requests` for HTTP calls and expects responses in
JSON format.
"""

from __future__ import annotations

from typing import Any, Iterable

import logging
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Configure module level logger
logger = logging.getLogger(__name__)

# Configure a session with retry and backoff for all HTTP requests
_retry = Retry(
    total=3,
    backoff_factor=1.0,
    status_forcelist=[500, 502, 503, 504],
    allowed_methods=["GET"],
)
_session = requests.Session()
_adapter = HTTPAdapter(max_retries=_retry)
_session.mount("http://", _adapter)
_session.mount("https://", _adapter)

# ----------------------------
# ChEMBL Target utilities
# ----------------------------

TARGET_FIELDS = [
    "pref_name",
    "target_chembl_id",
    "component_description",
    "component_id",
    "relationship",
    "gene",
    "uniprot_id",
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


def _parse_uniprot_id(xrefs: list[dict[str, str]]) -> str:
    """Extract a UniProt accession from a list of cross references.

    The ChEMBL API may label UniProt references with slightly different
    source database names; the check therefore normalises the label to
    upper case and matches common variants.
    """

    for x in xrefs:
        src = (x.get("xref_src_db") or "").upper()
        if src in {"UNIPROT", "UNIPROT ACCESSION", "UNIPROT ACC", "UNIPROTKB"}:
            ident = x.get("xref_id", "")
            if ident:
                return ident
    return ""


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


def _chunked(items: list[str], size: int) -> Iterable[list[str]]:
    """Yield successive ``size``-length chunks from ``items``.

    Parameters
    ----------
    items:
        Sequence to split into chunks.
    size:
        Maximum number of elements per chunk.

    Yields
    ------
    list[str]
        Slices of ``items`` with length up to ``size``.
    """
    if size <= 0:
        raise ValueError("size must be a positive integer")

    for i in range(0, len(items), size):
        yield items[i : i + size]


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

    gene = _parse_gene_synonyms(synonyms)
    ec_code = _parse_ec_codes(synonyms)
    alt_name = _parse_alt_names(synonyms)
    uniprot_id = _parse_uniprot_id(xrefs)
    hgnc_name, hgnc_id = _parse_hgnc(xrefs)

    return {
        "pref_name": data.get("pref_name", ""),
        "target_chembl_id": data.get("target_chembl_id", ""),
        "component_description": comp.get("component_description", ""),
        "component_id": comp.get("component_id", ""),
        "relationship": comp.get("relationship", ""),
        "gene": gene,
        "uniprot_id": uniprot_id,
        "chembl_alternative_name": alt_name,
        "ec_code": ec_code,
        "hgnc_name": hgnc_name,
        "hgnc_id": hgnc_id,
    }


def get_target(chembl_target_id: str) -> dict[str, Any]:
    """Fetch target data from ChEMBL for a single identifier.

    Parameters
    ----------
    chembl_target_id:
        ChEMBL target identifier.

    Returns
    -------
    dict
        A dictionary containing information about the target, including
        a ``uniprot_id`` when a UniProt cross reference is present. If the
        request fails an empty dictionary with pre-defined keys is returned.
    """
    if chembl_target_id in {"", "#N/A"}:
        return dict(EMPTY_TARGET)

    url = f"https://www.ebi.ac.uk/chembl/api/data/target/{chembl_target_id}?format=json"
    try:
        response = _session.get(url, timeout=30)
        response.raise_for_status()
    except requests.RequestException as exc:  # pragma: no cover - network
        logger.warning("Target request failed for %s: %s", chembl_target_id, exc)
        return dict(EMPTY_TARGET)

    try:
        data = response.json()
    except ValueError as exc:  # pragma: no cover - malformed JSON
        logger.warning("Failed to decode JSON for target %s: %s", chembl_target_id, exc)
        return dict(EMPTY_TARGET)
    return _parse_target_record(data)


def get_targets(ids: Iterable[str], chunk_size: int = 5) -> pd.DataFrame:
    """Fetch target records for ``ids``.

    Parameters
    ----------
    ids:
        Target identifiers to retrieve.
    chunk_size:
        Maximum number of IDs per HTTP request. ChEMBL rejects requests with
        very long query strings, so large input lists are split into chunks.

    Returns
    -------
    pandas.DataFrame
        DataFrame with one row per identifier. Failed chunks are skipped and
        logged; an entirely failed retrieval yields an empty DataFrame.
    """
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=TARGET_FIELDS)

    records: list[dict[str, Any]] = []
    for chunk in _chunked(valid, chunk_size):
        url = (
            "https://www.ebi.ac.uk/chembl/api/data/target.json?format=json&target_chembl_id__in="
            + ",".join(chunk)
        )
        try:
            response = _session.get(url, timeout=30)
            response.raise_for_status()
        except requests.RequestException as exc:  # pragma: no cover - network
            logger.warning("Bulk target request failed for %s: %s", chunk, exc)
            continue

        try:
            data = response.json()
        except ValueError as exc:  # pragma: no cover - malformed JSON
            logger.warning("Failed to decode JSON for targets %s: %s", chunk, exc)
            continue

        items = data.get("targets") or data.get("target") or []
        records.extend(_parse_target_record(item) for item in items)

    if not records:
        return pd.DataFrame(columns=TARGET_FIELDS)

    df = pd.DataFrame(records)
    return df.reindex(columns=TARGET_FIELDS)


# ----------------------------
# Assay utilities
# ----------------------------

ASSAY_URL = "https://www.ebi.ac.uk/chembl/api/data/assay/{id}?format=json"

ASSAY_COLUMNS = [
    "aidx",
    "assay_category",
    "assay_cell_type",
    "assay_chembl_id",
    "assay_classifications",
    "assay_group",
    "assay_organism",
    "assay_parameters",
    "assay_strain",
    "assay_subcellular_fraction",
    "assay_tax_id",
    "assay_test_type",
    "assay_tissue",
    "assay_type",
    "assay_type_description",
    "bao_format",
    "bao_label",
    "cell_chembl_id",
    "confidence_score",
    "description",
    "document_chembl_id",
    "src_assay_id",
    "src_id",
    "relationship_type",
    "target_chembl_id",
    "tissue_chembl_id",
    "variant_sequence.isoform",
    "variant_sequence.mutation",
    "variant_sequence.sequence",
]


def get_assay(chembl_assay_id: str) -> pd.DataFrame:
    """Retrieve assay information as a DataFrame.

    Parameters
    ----------
    chembl_assay_id:
        ChEMBL assay identifier.

    Returns
    -------
    pandas.DataFrame
        A single-row DataFrame with assay information.  If the request
        fails an empty DataFrame with predefined columns is returned.
    """
    if chembl_assay_id in {"", "#N/A"}:
        return pd.DataFrame(columns=ASSAY_COLUMNS)

    url = ASSAY_URL.format(id=chembl_assay_id)
    try:
        response = _session.get(url, timeout=30)
        response.raise_for_status()
    except requests.RequestException as exc:  # pragma: no cover - network
        logger.warning("Assay request failed for %s: %s", chembl_assay_id, exc)
        return pd.DataFrame(columns=ASSAY_COLUMNS)

    try:
        data = response.json()
    except ValueError as exc:  # pragma: no cover - malformed JSON
        logger.warning("Failed to decode JSON for assay %s: %s", chembl_assay_id, exc)
        return pd.DataFrame(columns=ASSAY_COLUMNS)
    df = pd.json_normalize(data)
    df = df.reindex(columns=ASSAY_COLUMNS)
    return df


def get_assays_all(ids: Iterable[str], chunk_size: int = 5) -> pd.DataFrame:
    """Fetch assay records for ``ids``.

    Parameters
    ----------
    ids:
        Assay identifiers to retrieve.
    chunk_size:
        Maximum number of IDs per HTTP request.

    Returns
    -------
    pandas.DataFrame
        Combined assay records.
    """
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=ASSAY_COLUMNS)

    records: list[pd.DataFrame] = []
    for chunk in _chunked(valid, chunk_size):
        url = (
            "https://www.ebi.ac.uk/chembl/api/data/assay.json?format=json&assay_chembl_id__in="
            + ",".join(chunk)
        )
        try:
            response = _session.get(url, timeout=1000)
            response.raise_for_status()
        except requests.RequestException as exc:  # pragma: no cover - network
            logger.warning("Bulk assay request failed for %s: %s", chunk, exc)
            continue

        try:
            data = response.json()
        except ValueError as exc:  # pragma: no cover - malformed JSON
            logger.warning("Failed to decode JSON for assays %s: %s", chunk, exc)
            continue

        items = data.get("assays") or data.get("assay") or []
        if items:
            # Normalise JSON then drop columns consisting solely of NaN values
            df_chunk = pd.json_normalize(items).dropna(axis="columns", how="all")
            if not df_chunk.empty:
                records.append(df_chunk)

    # If every request failed or returned only empty records, yield an empty frame

    if not records:
        return pd.DataFrame(columns=ASSAY_COLUMNS)

    df = pd.concat(records, ignore_index=True)
    return df.reindex(columns=ASSAY_COLUMNS)

def get_assays_notNull(ids: Iterable[str], chunk_size: int = 5) -> pd.DataFrame:
    """Fetch assay records for ``ids``.

    Parameters
    ----------
    ids:
        Assay identifiers to retrieve.
    chunk_size:
        Maximum number of IDs per HTTP request.

    Returns
    -------
    pandas.DataFrame
        Combined assay records.
    """
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=ASSAY_COLUMNS)

    records: list[pd.DataFrame] = []
    for chunk in _chunked(valid, chunk_size):
        url = (
            "https://www.ebi.ac.uk/chembl/api/data/assay.json?format=json&variant_sequence__isnull=false&assay_chembl_id__in="
            + ",".join(chunk)
        )
        try:
            response = _session.get(url, timeout=1000)
            response.raise_for_status()
        except requests.RequestException as exc:  # pragma: no cover - network
            logger.warning("Bulk assay request failed for %s: %s", chunk, exc)
            continue

        try:
            data = response.json()
        except ValueError as exc:  # pragma: no cover - malformed JSON
            logger.warning("Failed to decode JSON for assays %s: %s", chunk, exc)
            continue

        items = data.get("assays") or data.get("assay") or []
        if items:
            # Normalise JSON then drop columns consisting solely of NaN values
            df_chunk = pd.json_normalize(items).dropna(axis="columns", how="all")
            if not df_chunk.empty:
                records.append(df_chunk)

    # If every request failed or returned only empty records, yield an empty frame

    if not records:
        return pd.DataFrame(columns=ASSAY_COLUMNS)

    df = pd.concat(records, ignore_index=True)
    return df.reindex(columns=ASSAY_COLUMNS)
# ----------------------------
# Activity utilities
# ----------------------------

ACTIVITY_URL = "https://www.ebi.ac.uk/chembl/api/data/activity/{id}?format=json"

ACTIVITY_COLUMNS = [
    "activity_id",
    "assay_chembl_id",
    "document_chembl_id",
    "molecule_chembl_id",
    "standard_type",
    "standard_value",
    "standard_units",
    "standard_relation",
    "pchembl_value",
    "activity_comment",
    "data_validity_comment",
    "potential_duplicate",
    "bao_label",
    "src_id",
    "src_assay_id",
]


def get_activities(
    ids: Iterable[str],
    chunk_size: int = 5,
    timeout: float = 30.0,
) -> pd.DataFrame:
    """Fetch activity records for ``ids``.

    Parameters
    ----------
    ids:
        Activity identifiers to retrieve.
    chunk_size:
        Maximum number of IDs per HTTP request.
    timeout:
        Timeout in seconds for each HTTP request.

    Returns
    -------
    pandas.DataFrame
        Combined activity records.
    """
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=ACTIVITY_COLUMNS)

    records: list[pd.DataFrame] = []
    for chunk in _chunked(valid, chunk_size):
        url = (
            "https://www.ebi.ac.uk/chembl/api/data/activity.json?format=json&activity_id__in="
            + ",".join(chunk)
        )
        try:
            response = _session.get(url, timeout=timeout)
            response.raise_for_status()
        except requests.RequestException as exc:  # pragma: no cover - network
            logger.warning("Bulk activity request failed for %s: %s", chunk, exc)
            continue

        try:
            data = response.json()
        except ValueError as exc:  # pragma: no cover - malformed JSON
            logger.warning("Failed to decode JSON for activities %s: %s", chunk, exc)
            continue

        items = data.get("activities") or data.get("activity") or []
        if items:
            records.append(pd.json_normalize(items))

    if not records:
        return pd.DataFrame(columns=ACTIVITY_COLUMNS)

    df = pd.concat(records, ignore_index=True)
    return df.reindex(columns=ACTIVITY_COLUMNS)


# ----------------------------
# Test item (compound) utilities
# ----------------------------

TESTITEM_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule/{id}?format=json"

TESTITEM_COLUMNS = [
    "molecule_chembl_id",
    "pref_name",
    "max_phase",
    "molecule_type",
    "first_approval",
    "oral",
    "parenteral",
    "topical",
    "black_box_warning",
    "structure_type",
    "molecule_structures.canonical_smiles",
    "molecule_structures.standard_inchi",
    "molecule_structures.standard_inchi_key",
]


def get_testitem(ids: Iterable[str], chunk_size: int = 5) -> pd.DataFrame:
    """Fetch compound records for ``ids``.

    Parameters
    ----------
    ids:
        Molecule identifiers to retrieve.
    chunk_size:
        Maximum number of IDs per HTTP request.

    Returns
    -------
    pandas.DataFrame
        Combined compound records.
    """
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=TESTITEM_COLUMNS)

    records: list[pd.DataFrame] = []
    for chunk in _chunked(valid, chunk_size):
        url = (
            "https://www.ebi.ac.uk/chembl/api/data/molecule.json?format=json&molecule_chembl_id__in="
            + ",".join(chunk)
        )
        try:
            response = _session.get(url, timeout=30)
            response.raise_for_status()
        except requests.RequestException as exc:  # pragma: no cover - network
            logger.warning("Bulk molecule request failed for %s: %s", chunk, exc)
            continue

        try:
            data = response.json()
        except ValueError as exc:  # pragma: no cover - malformed JSON
            logger.warning("Failed to decode JSON for molecules %s: %s", chunk, exc)
            continue

        items = data.get("molecules") or data.get("molecule") or []
        if items:
            records.append(pd.json_normalize(items))
 # Drop empty or all-NA frames to avoid deprecation warnings in pandas
    records = [r for r in records if not r.empty and not r.isna().all().all()]
    if not records:
        return pd.DataFrame(columns=TESTITEM_COLUMNS)

    df = pd.concat(records, ignore_index=True)
    return df.reindex(columns=TESTITEM_COLUMNS)


# ----------------------------
# Document utilities
# ----------------------------

DOCUMENT_URL = "https://www.ebi.ac.uk/chembl/api/data/document/{id}?format=json"

DOCUMENT_COLUMNS = [
    "document_chembl_id",
    "title",
    "abstract",
    "doi",
    "year",
    "journal",
    "journal_abbrev",
    "volume",
    "issue",
    "first_page",
    "last_page",
    "pubmed_id",
    "authors",
    "source",
]


def get_document(chembl_document_id: str) -> pd.DataFrame:
    """Retrieve document information for a single identifier.

    Parameters
    ----------
    chembl_document_id:
        ChEMBL document identifier.

    Returns
    -------
    pandas.DataFrame
        Single-row DataFrame with document data.  An empty DataFrame is
        returned if the request fails or the ID is blank.
    """
    if chembl_document_id in {"", "#N/A"}:
        return pd.DataFrame(columns=DOCUMENT_COLUMNS)

    url = DOCUMENT_URL.format(id=chembl_document_id)
    try:
        response = _session.get(url, timeout=30)
        response.raise_for_status()
    except requests.RequestException as exc:  # pragma: no cover - network
        logger.warning("Document request failed for %s: %s", chembl_document_id, exc)
        return pd.DataFrame(columns=DOCUMENT_COLUMNS)

    try:
        data = response.json()
    except ValueError as exc:  # pragma: no cover - malformed JSON
        logger.warning(
            "Failed to decode JSON for document %s: %s", chembl_document_id, exc
        )
        return pd.DataFrame(columns=DOCUMENT_COLUMNS)

    df = pd.json_normalize(data)
    return df.reindex(columns=DOCUMENT_COLUMNS)


def get_documents(ids: Iterable[str], chunk_size: int = 5) -> pd.DataFrame:
    """Fetch document records for ``ids``.

    Parameters
    ----------
    ids:
        Document identifiers to retrieve.
    chunk_size:
        Maximum number of IDs per HTTP request.

    Returns
    -------
    pandas.DataFrame
        Combined document records.
    """
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=DOCUMENT_COLUMNS)

    records: list[pd.DataFrame] = []
    for chunk in _chunked(valid, chunk_size):
        url = (
            "https://www.ebi.ac.uk/chembl/api/data/document.json?format=json&document_chembl_id__in="
            + ",".join(chunk)
        )
        try:
            response = _session.get(url, timeout=30)
            response.raise_for_status()
        except requests.RequestException as exc:  # pragma: no cover - network
            logger.warning("Bulk document request failed for %s: %s", chunk, exc)
            continue

        try:
            data = response.json()
        except ValueError as exc:  # pragma: no cover - malformed JSON
            logger.warning("Failed to decode JSON for documents %s: %s", chunk, exc)
            continue

        items = data.get("documents") or data.get("document") or []
        if items:
            records.append(pd.json_normalize(items))

    if not records:
        return pd.DataFrame(columns=DOCUMENT_COLUMNS)

    df = pd.concat(records, ignore_index=True)
    return df.reindex(columns=DOCUMENT_COLUMNS)


def extend_target(
    df: pd.DataFrame, chembl_column: str = "task_chembl_id", chunk_size: int = 5
) -> pd.DataFrame:
    """Augment a DataFrame with target information.

    Parameters
    ----------
    df:
        Input table containing a column with ChEMBL target IDs.
    chembl_column:
        Name of the column holding the target IDs. Default is
        ``"task_chembl_id"``.
    chunk_size:
        Maximum number of IDs per HTTP request when fetching targets.

    Returns
    -------
    pandas.DataFrame
        Original data combined with the expanded target information.
    """
    if chembl_column not in df.columns:
        raise ValueError(f"column '{chembl_column}' not found in DataFrame")

    ids = df[chembl_column].astype(str).tolist()
    targets = get_targets(ids, chunk_size=chunk_size)
    merged = df.merge(
        targets,
        how="left",
        left_on=chembl_column,
        right_on="target_chembl_id",
    )
    return merged.rename(
        columns={
            "pref_name": "chembl_pref_name",
            "target_chembl_id": "chembl_target_chembl_id",
            "component_description": "chembl_component_description",
            "component_id": "chembl_component_id",
            "relationship": "chembl_relationship",
            "gene": "chembl_gene",
            "chembl_alternative_name": "chembl_alternative_name",
            "ec_code": "chembl_ec_code",
            "hgnc_name": "chembl_hgnc_name",
            "hgnc_id": "chembl_hgnc_id",
        }
    )