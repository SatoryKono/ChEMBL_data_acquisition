"""Assay and activity helpers for the ChEMBL API."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from urllib.parse import urljoin

import pandas as pd

from .clients.chembl_client import ChemblClient, _chunked
from .config import ApiCfg
from .log import logger
from .io.pandas_utils import json_normalize_pyarrow

ASSAY_VARIANT_COLUMN_ALIASES = {
    "variant_sequence.isoform": "isoform",
    "variant_sequence.mutation": "mutation",
    "variant_sequence.sequence": "sequence",
}


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
    "isoform",
    "mutation",
    "sequence",
]


def _apply_assay_column_aliases(df: pd.DataFrame) -> pd.DataFrame:
    """Rename nested variant sequence columns to flattened aliases."""
    if df.empty:
        return df
    return df.rename(columns=ASSAY_VARIANT_COLUMN_ALIASES)


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
    "standard_lower_value",
    "standard_value",
    "standard_upper_value",
    "standard_units",
    "standard_relation",
    "type",
    "value",
    "upper_value",
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

TESTITEM_STRUCTURE_COLUMNS = {
    "molecule_structures.canonical_smiles": "canonical_smiles",
    "molecule_structures.standard_inchi": "standard_inchi",
    "molecule_structures.standard_inchi_key": "standard_inchi_key",
}

TESTITEM_COLUMNS = [
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
    "pref_name",
    "max_phase",
    "molecule_type",
    "first_approval",
    "oral",
    "parenteral",
    "topical",
    "black_box_warning",
    "structure_type",
    "canonical_smiles",
    "standard_inchi",
    "standard_inchi_key",
]


def get_assay(
    chembl_assay_id: str,
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Retrieve assay information as a DataFrame.

    Parameters
    ----------
    chembl_assay_id:
        Identifier of the assay to fetch.
    cfg:
        API configuration providing base URL and timeouts.
    client:
        ChemblClient used for HTTP requests and caching.
    timeout:
        Optional override for the read timeout in seconds.

    Returns
    -------
    pandas.DataFrame
        Normalised assay information or an empty frame when ``chembl_assay_id``
        is empty or the record is missing.
    """
    if chembl_assay_id in {"", "#N/A"}:
        return pd.DataFrame(columns=ASSAY_COLUMNS)
    base = cfg.chembl_base.rstrip("/")
    url = f"{base}/assay/{chembl_assay_id}?format=json"
    effective_timeout = timeout if timeout is not None else cfg.timeout_read
    data = client.request_json(url, cfg=cfg, timeout=effective_timeout)
    items = data.get("assays") or data.get("assay") or []
    if not items:
        return pd.DataFrame(columns=ASSAY_COLUMNS)
    df = json_normalize_pyarrow(items).dropna(axis="columns", how="all")
    df = _apply_assay_column_aliases(df)
    return df.reindex(columns=ASSAY_COLUMNS)


def get_assays(
    ids: Iterable[str],
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    chunk_size: int = 5,
    timeout: float | None = None,
    require_variant_sequence: bool = False,
) -> pd.DataFrame:
    """Fetch assay records for *ids*.

    Parameters
    ----------
    ids:
        Assay identifiers to retrieve.
    cfg:
        API configuration providing base URL and timeouts.
    client:
        ChemblClient used for HTTP requests and caching.
    chunk_size:
        Maximum number of IDs per HTTP request.
    timeout:
        Optional override for the read timeout in seconds.
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
    base = f"{cfg.chembl_base.rstrip('/')}/assay.json?format=json"
    if require_variant_sequence:
        base += "&variant_sequence__isnull=false"
    effective_timeout = timeout if timeout is not None else cfg.timeout_read
    for chunk in _chunked(valid, chunk_size):
        chunk_key = ",".join(chunk)
        logger.info(
            "chunk_start", extra={"stage": "chunk_start", "chunk_key": chunk_key}
        )
        url = f"{base}&assay_chembl_id__in={chunk_key}&limit={len(chunk)}"
        chunk_frames: list[pd.DataFrame] = []
        next_url: str | None = url
        while next_url:
            data = client.request_json(next_url, cfg=cfg, timeout=effective_timeout)
            items = data.get("assays") or data.get("assay") or []
            if items:
                df_chunk = json_normalize_pyarrow(items).dropna(
                    axis="columns", how="all"
                )
                if not df_chunk.empty:
                    chunk_frames.append(_apply_assay_column_aliases(df_chunk))
            page_meta = data.get("page_meta") or {}
            next_token = page_meta.get("next")
            next_url = urljoin(cfg.chembl_base, next_token) if next_token else None
        if chunk_frames:
            records.append(pd.concat(chunk_frames, ignore_index=True))
            logger.info(
                "chunk_done",
                extra={"stage": "chunk_done", "chunk_key": chunk_key},
            )
            continue
        logger.info("chunk_skip", extra={"stage": "chunk_skip", "chunk_key": chunk_key})
    if not records:
        return pd.DataFrame(columns=ASSAY_COLUMNS)
    df = pd.concat(records, ignore_index=True)
    df = _apply_assay_column_aliases(df)
    return df.reindex(columns=ASSAY_COLUMNS)


def get_activities(
    ids: Iterable[str],
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    chunk_size: int = 5,
    timeout: float | None = None,
    extra_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Fetch activity records for *ids*.

    Parameters
    ----------
    ids:
        Activity identifiers to retrieve.
    cfg:
        API configuration providing base URL and timeouts.
    client:
        ChemblClient used for HTTP requests and caching.
    chunk_size:
        Maximum number of IDs per HTTP request.
    timeout:
        Optional override for the read timeout in seconds.

    Returns
    -------
    pandas.DataFrame
        Combined activity records.
    """
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=ACTIVITY_COLUMNS)

    records: list[pd.DataFrame] = []
    base = f"{cfg.chembl_base.rstrip('/')}/activity.json?format=json"
    effective_timeout = timeout if timeout is not None else cfg.timeout_read
    for chunk in _chunked(valid, chunk_size):
        chunk_key = ",".join(chunk)
        logger.info(
            "chunk_start", extra={"stage": "chunk_start", "chunk_key": chunk_key}
        )
        url = f"{base}&activity_id__in={chunk_key}"
        data = client.request_json(url, cfg=cfg, timeout=effective_timeout)
        items = data.get("activities") or data.get("activity") or []
        if items:
            records.append(json_normalize_pyarrow(items))
            logger.info(
                "chunk_done", extra={"stage": "chunk_done", "chunk_key": chunk_key}
            )
        else:
            logger.info(
                "chunk_skip", extra={"stage": "chunk_skip", "chunk_key": chunk_key}
            )
    if not records:
        columns = list(ACTIVITY_COLUMNS)
        if extra_columns:
            for column in extra_columns:
                if column not in columns:
                    columns.append(column)
        return pd.DataFrame(columns=columns)
    df = pd.concat(records, ignore_index=True)
    columns = list(ACTIVITY_COLUMNS)
    if extra_columns:
        for column in extra_columns:
            if column not in columns:
                columns.append(column)
    return df.reindex(columns=columns)


def get_testitem(
    ids: Iterable[str],
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    chunk_size: int = 5,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Fetch compound records for *ids*.

    Parameters
    ----------
    ids:
        Molecule identifiers to retrieve.
    cfg:
        API configuration providing base URL and timeouts.
    client:
        ChemblClient used for HTTP requests and caching.
    chunk_size:
        Maximum number of IDs per HTTP request.
    timeout:
        Optional override for the read timeout in seconds.

    Returns
    -------
    pandas.DataFrame
        Combined compound records.
    """
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=TESTITEM_COLUMNS)

    records: list[pd.DataFrame] = []
    base = f"{cfg.chembl_base.rstrip('/')}/molecule.json?format=json"
    effective_timeout = timeout if timeout is not None else cfg.timeout_read
    for chunk in _chunked(valid, chunk_size):
        chunk_key = ",".join(chunk)
        logger.info(
            "chunk_start", extra={"stage": "chunk_start", "chunk_key": chunk_key}
        )
        url = f"{base}&molecule_chembl_id__in={chunk_key}"
        data = client.request_json(url, cfg=cfg, timeout=effective_timeout)
        items = data.get("molecules") or data.get("molecule") or []
        if items:
            records.append(json_normalize_pyarrow(items))
            logger.info(
                "chunk_done", extra={"stage": "chunk_done", "chunk_key": chunk_key}
            )
        else:
            logger.info(
                "chunk_skip", extra={"stage": "chunk_skip", "chunk_key": chunk_key}
            )
    if not records:
        return pd.DataFrame(columns=TESTITEM_COLUMNS)
    df = pd.concat(records, ignore_index=True)
    df = df.rename(columns=TESTITEM_STRUCTURE_COLUMNS)
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
