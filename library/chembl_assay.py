"""Assay and activity helpers for the ChEMBL API."""

from __future__ import annotations

from typing import Iterable

import logging

import pandas as pd

from .chembl_client import _chunked, request_json
from .config import Config

logger = logging.getLogger(__name__)

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


ACTIVITY_COLUMNS = [
    "activity_id",
    "assay_chembl_id",
    "document_chembl_id",
    "molecule_chembl_id",
    "record_id",
    "assay_description",
    "bao_label",
    "bao_format",
    "standard_type",
    "standard_value",
    "standard_units",
    "standard_relation",
    "type",
    "value",
    "units",
    "relation",
    "qudt_units",
    "pchembl_value",
    "activity_comment",
    "data_validity_comment",
    "data_validity_description",
    "potential_duplicate",
    "src_id",
    "src_assay_id",
    "assay_variant_accession",
    "assay_variant_mutation",
]

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


def get_assay(
    chembl_assay_id: str,
    *,
    cfg: Config,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Retrieve assay information as a DataFrame."""
    if chembl_assay_id in {"", "#N/A"}:
        return pd.DataFrame(columns=ASSAY_COLUMNS)
    url = ASSAY_URL.format(id=chembl_assay_id)
    data = request_json(cfg, url, timeout=timeout)
    items = data.get("assays") or data.get("assay") or []
    if not items:
        return pd.DataFrame(columns=ASSAY_COLUMNS)
    df = pd.json_normalize(items, dtype_backend="pyarrow").dropna(  # type: ignore[call-arg]
        axis="columns", how="all"
    )
    return df.reindex(columns=ASSAY_COLUMNS)


def get_assays(
    ids: Iterable[str],
    *,
    cfg: Config,
    chunk_size: int = 5,
    timeout: float | None = None,
    require_variant_sequence: bool = False,
) -> pd.DataFrame:
    """Fetch assay records for *ids*.

    Parameters
    ----------
    ids:
        Assay identifiers to retrieve.
    chunk_size:
        Maximum number of IDs per HTTP request.
    timeout:
        Timeout in seconds for each HTTP request.
    require_variant_sequence:
        If ``True``, only assays with a non-null ``variant_sequence`` are returned.

    Returns
    -------
    pandas.DataFrame
        Combined assay records.

    """
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=ASSAY_COLUMNS)

    records: list[pd.DataFrame] = []
    base = "https://www.ebi.ac.uk/chembl/api/data/assay.json?format=json"
    if require_variant_sequence:
        base += "&variant_sequence__isnull=false"
    for chunk in _chunked(valid, chunk_size):
        url = f"{base}&assay_chembl_id__in={','.join(chunk)}"
        data = request_json(cfg, url, timeout=timeout)
        items = data.get("assays") or data.get("assay") or []
        if items:
            df_chunk = pd.json_normalize(items, dtype_backend="pyarrow").dropna(  # type: ignore[call-arg]
                axis="columns", how="all"
            )
            if not df_chunk.empty:
                records.append(df_chunk)
    if not records:
        return pd.DataFrame(columns=ASSAY_COLUMNS)
    df = pd.concat(records, ignore_index=True)
    return df.reindex(columns=ASSAY_COLUMNS)


def get_activities(
    ids: Iterable[str],
    *,
    cfg: Config,
    chunk_size: int = 5,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Fetch activity records for *ids*.

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
    base = "https://www.ebi.ac.uk/chembl/api/data/activity.json?format=json"
    for chunk in _chunked(valid, chunk_size):
        url = f"{base}&activity_id__in={','.join(chunk)}"
        data = request_json(cfg, url, timeout=timeout)
        items = data.get("activities") or data.get("activity") or []
        if items:
            records.append(pd.json_normalize(items, dtype_backend="pyarrow"))  # type: ignore[call-arg]
    if not records:
        return pd.DataFrame(columns=ACTIVITY_COLUMNS)
    df = pd.concat(records, ignore_index=True)
    return df.reindex(columns=ACTIVITY_COLUMNS)


def get_testitem(
    ids: Iterable[str],
    *,
    cfg: Config,
    chunk_size: int = 5,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Fetch compound records for *ids*."""
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=TESTITEM_COLUMNS)

    records: list[pd.DataFrame] = []
    base = "https://www.ebi.ac.uk/chembl/api/data/molecule.json?format=json"
    for chunk in _chunked(valid, chunk_size):
        url = f"{base}&molecule_chembl_id__in={','.join(chunk)}"
        data = request_json(cfg, url, timeout=timeout)
        items = data.get("molecules") or data.get("molecule") or []
        if items:
            records.append(pd.json_normalize(items, dtype_backend="pyarrow"))  # type: ignore[call-arg]
    if not records:
        return pd.DataFrame(columns=TESTITEM_COLUMNS)
    df = pd.concat(records, ignore_index=True)
    return df.reindex(columns=TESTITEM_COLUMNS)


__all__ = [
    "get_assay",
    "get_assays",
    "get_activities",
    "get_testitem",
    "ASSAY_COLUMNS",
    "ACTIVITY_COLUMNS",
    "TESTITEM_COLUMNS",
]
