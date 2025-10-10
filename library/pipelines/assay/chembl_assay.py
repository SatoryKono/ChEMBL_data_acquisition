"""Assay and activity helpers for the ChEMBL API."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from urllib.parse import urlencode, urljoin

import pandas as pd
import requests

from library.clients import ChemblClient, _chunked

from ...common.log import logger
from ...common.pandas_utils import json_normalize_pyarrow
from ...config import TESTITEM_FIELD_DEFAULTS, ApiCfg

# ChEMBL bulk endpoints accept at most 25 identifiers per request when using
# ``__in`` filters. Requests above this threshold return HTTP 414 "URI Too
# Long", so keep the chunk size at or below this value to avoid retries.
MAX_ASSAY_CHUNK_SIZE = 25

# ``_ASSAY_MAX_IDS_PER_REQUEST`` used to be the public name for the chunk size
# limit.  Keep the alias for compatibility with older call sites that still
# reference it indirectly (for example through ``eval`` of serialized
# expressions).  The leading underscore mirrors the historical constant name and
# avoids polluting the public namespace beyond backward compatibility needs.
_ASSAY_MAX_IDS_PER_REQUEST = MAX_ASSAY_CHUNK_SIZE

MAX_ACTIVITY_CHUNK_SIZE = 20

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
    "activity_chembl_id",
    "assay_chembl_id",
    "target_chembl_id",
    "document_chembl_id",
    "molecule_chembl_id",
    "salt_chembl_id",
    "record_id",
    "assay_description",
    "assay_variant_accession",
    "assay_variant_mutation",
    "bao_label",
    "bao_format",
    "bao_endpoint",
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
    "log_value",
    "activity_comment",
    "data_validity_comment",
    "data_validity_description",
    "activity_properties",
    "properties_hash",
    "potential_duplicate",
    "pipeline_version",
    "timestamp_utc",
    "compound_key",
    "compound_name",
    "multmol_assay",
    "approx_cited_activity",
    "shuffled_cit",
    "exact_cited_activity",
    "higly_correlated_cit",
    "review_doc",
    "rounded_data_citation",
    "original_activity_approx",
    "original_activity_exact",
    "nstereo",
    "lower_value",
    "src_id",
    "src_assay_id",
    "action_type",
]


_ACTIVITY_EXTRA_FIELDS: tuple[str, ...] = ("standard_text_value",)

ACTIVITY_QUERY_FIELDS: tuple[str, ...] = tuple(
    dict.fromkeys(list(ACTIVITY_COLUMNS) + list(_ACTIVITY_EXTRA_FIELDS))
)

TESTITEM_PUBCHEM_COLUMNS: tuple[str, ...] = (
    "pubchem_cid",
    "pubchem_iupac_name",
    "pubchem_molecular_formula",
    "pubchem_isomeric_smiles",
    "pubchem_canonical_smiles",
    "pubchem_inchi",
    "pubchem_inchikey",
)


TESTITEM_STRUCTURE_COLUMNS = {
    "pubchem.cid": "pubchem_cid",
    "pubchem.canonical_smiles": "pubchem_canonical_smiles",
    "pubchem.isomeric_smiles": "pubchem_isomeric_smiles",
    "pubchem.inchi": "pubchem_inchi",
    "pubchem.inchikey": "pubchem_inchikey",
    "pubchem.inchi_key": "pubchem_inchikey",
    "pubchem.iupac_name": "pubchem_iupac_name",
    "pubchem.molecular_formula": "pubchem_molecular_formula",
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
    *TESTITEM_PUBCHEM_COLUMNS,
]


TESTITEM_QUERY_FIELDS = TESTITEM_FIELD_DEFAULTS


# ChEMBL API rejects URLs above ~4096 characters, keep a conservative margin.
MAX_TESTITEM_URL_LENGTH = 4000


def _split_chunk_for_url(
    chunk: Sequence[str],
    base_url: str,
    *,
    max_length: int = MAX_TESTITEM_URL_LENGTH,
) -> Iterator[list[str]]:
    """Yield sub-chunks that keep the request URL within ``max_length``."""

    prefix = f"{base_url}&molecule_chembl_id__in="
    prefix_length = len(prefix)
    if prefix_length >= max_length:
        raise ValueError("base URL exceeds maximum request length")

    buffer: list[str] = []
    buffer_length = 0
    for identifier in chunk:
        separator = 1 if buffer else 0
        candidate_length = buffer_length + separator + len(identifier)
        if prefix_length + candidate_length > max_length and buffer:
            yield buffer
            buffer = [identifier]
            buffer_length = len(identifier)
        elif prefix_length + candidate_length > max_length:
            yield [identifier]
            buffer = []
            buffer_length = 0
        else:
            buffer.append(identifier)
            buffer_length = candidate_length

    if buffer:
        yield buffer


def _fetch_testitem_chunk(
    identifiers: Sequence[str],
    *,
    base: str,
    cfg: ApiCfg,
    client: ChemblClient,
    timeout: float,
) -> list[pd.DataFrame]:
    """Fetch data for ``identifiers`` with adaptive splitting on HTTP 400/timeouts."""

    pending: list[list[str]] = [list(identifiers)]
    frames: list[pd.DataFrame] = []
    while pending:
        current = pending.pop()
        if not current:
            continue
        chunk_key = ",".join(current)
        logger.info(
            "chunk_start", extra={"stage": "chunk_start", "chunk_key": chunk_key}
        )
        url = f"{base}&molecule_chembl_id__in={chunk_key}"
        next_url: str | None = url
        seen_urls: set[str] = set()
        chunk_frames: list[pd.DataFrame] = []
        try:
            while next_url:
                absolute_url = urljoin(cfg.chembl_base, next_url)
                if absolute_url in seen_urls:
                    logger.warning(
                        "pagination_loop_detected",
                        extra={"stage": "chunk_loop", "chunk_key": chunk_key},
                    )
                    break
                seen_urls.add(absolute_url)
                data = client.request_json(absolute_url, cfg=cfg, timeout=timeout)
                items = data.get("molecules") or data.get("molecule") or []
                if items:
                    chunk_frames.append(json_normalize_pyarrow(items))
                page_meta = data.get("page_meta") or {}
                next_token = page_meta.get("next")
                next_url = urljoin(cfg.chembl_base, next_token) if next_token else None
        except requests.RequestException as exc:
            response = getattr(exc, "response", None)
            status_code = response.status_code if response is not None else None
            if status_code == 400 and len(current) > 1:
                midpoint = max(1, len(current) // 2)
                first = current[:midpoint]
                second = current[midpoint:]
                logger.warning(
                    "chunk_http_400_split",
                    extra={
                        "stage": "chunk_split",
                        "chunk_key": chunk_key,
                        "size": len(current),
                        "status": status_code,
                    },
                )
                if second:
                    pending.append(second)
                if first:
                    pending.append(first)
                continue
            if _is_retryable_chunk_error(exc) and len(current) > 1:
                midpoint = max(1, len(current) // 2)
                first = current[:midpoint]
                second = current[midpoint:]
                logger.warning(
                    "chunk_retry_split",
                    extra={
                        "stage": "chunk_split",
                        "chunk_key": chunk_key,
                        "size": len(current),
                        "error": exc.__class__.__name__,
                    },
                )
                if second:
                    pending.append(second)
                if first:
                    pending.append(first)
                continue
            raise
        if chunk_frames:
            frames.append(pd.concat(chunk_frames, ignore_index=True))
            logger.info(
                "chunk_done", extra={"stage": "chunk_done", "chunk_key": chunk_key}
            )
        else:
            logger.info(
                "chunk_skip", extra={"stage": "chunk_skip", "chunk_key": chunk_key}
            )
    return frames


def _is_retryable_chunk_error(exc: requests.RequestException) -> bool:
    """Return ``True`` when *exc* indicates a transient connection failure."""

    if isinstance(exc, requests.Timeout):
        return True
    message_parts = [str(exc)]
    cause = getattr(exc, "__cause__", None)
    context = getattr(exc, "__context__", None)
    for detail in (cause, context):
        if detail is not None:
            message_parts.append(str(detail))
    combined = " ".join(part for part in message_parts if part)
    lowered = combined.lower()
    return any(
        token in lowered
        for token in ("timed out", "timeout", "temporarily unavailable")
    )


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
    base_root = cfg.chembl_base.rstrip("/")
    base_url = f"{base_root}/assay.json"
    base_query: dict[str, str] = {"format": "json"}
    if require_variant_sequence:
        base_query["variant_sequence__isnull"] = "false"
    join_base = f"{base_root}/"

    def _build_url(extra: dict[str, object]) -> str:
        params = base_query | {k: str(v) for k, v in extra.items()}
        return f"{base_url}?{urlencode(params)}"

    def _fetch_chunk(url: str) -> list[pd.DataFrame]:
        frames: list[pd.DataFrame] = []
        next_url: str | None = url
        while next_url:
            data = client.request_json(next_url, cfg=cfg, timeout=effective_timeout)
            items = data.get("assays") or data.get("assay") or []
            if items:
                df_chunk = json_normalize_pyarrow(items).dropna(
                    axis="columns", how="all"
                )
                if not df_chunk.empty:
                    frames.append(_apply_assay_column_aliases(df_chunk))
            page_meta = data.get("page_meta") or {}
            next_token = page_meta.get("next")
            next_url = urljoin(join_base, next_token) if next_token else None
        return frames

    def _fallback_single(identifier: str) -> list[pd.DataFrame]:
        fallback_df = get_assay(
            identifier,
            cfg=cfg,
            client=client,
            timeout=effective_timeout,
        )
        if require_variant_sequence and not fallback_df.empty:
            sequence_series = fallback_df.get("sequence")
            if sequence_series is not None:
                fallback_df = fallback_df[sequence_series.notna()]
        if fallback_df.empty:
            return []
        return [fallback_df]

    def _filter_variant_frames(frames: list[pd.DataFrame]) -> list[pd.DataFrame]:
        if not require_variant_sequence:
            return frames
        filtered_frames: list[pd.DataFrame] = []
        for frame in frames:
            sequence_series = frame.get("sequence")
            if sequence_series is None:
                continue
            filtered = frame[sequence_series.notna()]
            if not filtered.empty:
                filtered_frames.append(filtered)
        return filtered_frames

    def _fetch_single(identifier: str) -> list[pd.DataFrame]:
        single_url = _build_url({"assay_chembl_id": identifier, "limit": 1})
        try:
            frames = _fetch_chunk(single_url)
        except requests.HTTPError as exc:  # pragma: no cover - exercised via caller
            response = exc.response
            if response is not None and response.status_code == 404:
                logger.info(
                    "assay_missing",
                    extra={"stage": "chunk_skip", "assay_chembl_id": identifier},
                )
                return _fallback_single(identifier)
            raise
        except requests.RequestException as exc:
            logger.warning(
                "single_fetch_retry",
                extra={
                    "stage": "chunk_retry",
                    "assay_chembl_id": identifier,
                    "error": exc.__class__.__name__,
                },
            )
            return _fallback_single(identifier)
        else:
            frames = _filter_variant_frames(frames)
            return frames

    effective_timeout = timeout if timeout is not None else cfg.timeout_read
    effective_chunk_size = min(int(chunk_size), MAX_ASSAY_CHUNK_SIZE)
    if effective_chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    if effective_chunk_size < chunk_size:
        logger.debug(
            "assay_chunk_clamped",
            extra={
                "requested_chunk_size": chunk_size,
                "effective_chunk_size": effective_chunk_size,
                "limit": MAX_ASSAY_CHUNK_SIZE,
                "stage": "chunk_prepare",
            },
        )

    for chunk in _chunked(valid, effective_chunk_size):
        chunk_key = ",".join(chunk)
        logger.info(
            "chunk_start", extra={"stage": "chunk_start", "chunk_key": chunk_key}
        )
        url = _build_url({"assay_chembl_id__in": ",".join(chunk), "limit": len(chunk)})
        try:
            chunk_frames = _fetch_chunk(url)
        except requests.HTTPError as exc:
            response = exc.response
            if response is not None and response.status_code == 404:
                logger.warning(
                    "chunk_split_retry",
                    extra={
                        "stage": "chunk_retry",
                        "chunk_key": chunk_key,
                        "chunk_size": len(chunk),
                        "status": 404,
                    },
                )
                chunk_frames = []
                for identifier in chunk:
                    chunk_frames.extend(_fetch_single(identifier))
            else:
                raise
        except requests.RequestException as exc:
            logger.warning(
                "chunk_network_retry",
                extra={
                    "stage": "chunk_retry",
                    "chunk_key": chunk_key,
                    "chunk_size": len(chunk),
                    "error": exc.__class__.__name__,
                },
            )
            chunk_frames = []
            for identifier in chunk:
                try:
                    chunk_frames.extend(_fetch_single(identifier))
                except requests.RequestException:
                    raise
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
    effective_chunk_size = max(1, min(int(chunk_size), MAX_ACTIVITY_CHUNK_SIZE))
    if effective_chunk_size != chunk_size:
        logger.debug(
            "activity_chunk_size_clamped",
            extra={
                "requested": int(chunk_size),
                "effective": effective_chunk_size,
                "limit": MAX_ACTIVITY_CHUNK_SIZE,
            },
        )
    base_root = cfg.chembl_base.rstrip("/")
    base_url = f"{base_root}/activity.json"
    base_params: list[tuple[str, str]] = [("format", "json")]
    if ACTIVITY_QUERY_FIELDS:
        base_params.append(("fields", ",".join(ACTIVITY_QUERY_FIELDS)))
    join_base = f"{base_root}/"

    def _build_url(extra: dict[str, object]) -> str:
        params = base_params + [(key, str(value)) for key, value in extra.items()]
        return f"{base_url}?{urlencode(params)}"

    def _fetch_paginated(url: str) -> list[pd.DataFrame]:
        frames: list[pd.DataFrame] = []
        next_url: str | None = url
        while next_url:
            data = client.request_json(next_url, cfg=cfg, timeout=effective_timeout)
            items = data.get("activities") or data.get("activity") or []
            if items:
                df_chunk = json_normalize_pyarrow(items)
                if not df_chunk.empty:
                    frames.append(df_chunk)
            page_meta = data.get("page_meta") or {}
            next_token = page_meta.get("next")
            next_url = urljoin(join_base, next_token) if next_token else None
        return frames

    def _fetch_single(identifier: str) -> list[pd.DataFrame]:
        if identifier in {"", "#N/A"}:
            return []
        try:
            frames = _fetch_paginated(
                _build_url({"activity_id": identifier, "limit": 1})
            )
        except requests.HTTPError as exc:
            response = exc.response
            status = getattr(response, "status_code", None)
            if status == 404:
                logger.info(
                    "activity_missing",
                    extra={"stage": "chunk_skip", "activity_id": identifier},
                )
                return []
            raise
        except requests.RequestException as exc:
            logger.warning(
                "single_fetch_network_skip",
                extra={
                    "stage": "chunk_retry",
                    "activity_id": identifier,
                    "error": exc.__class__.__name__,
                },
            )
            return []
        return frames

    effective_timeout = timeout if timeout is not None else cfg.timeout_read
    for chunk in _chunked(valid, effective_chunk_size):
        chunk_key = ",".join(chunk)
        logger.info(
            "chunk_start",
            extra={
                "stage": "chunk_start",
                "chunk_key": chunk_key,
                "chunk_size": len(chunk),
            },
        )
        url = _build_url({"activity_id__in": chunk_key, "limit": len(chunk)})
        try:
            chunk_frames = _fetch_paginated(url)
        except requests.HTTPError as exc:
            response = exc.response
            status = getattr(response, "status_code", None)
            if status == 404:
                logger.warning(
                    "chunk_split_retry",
                    extra={
                        "stage": "chunk_retry",
                        "chunk_key": chunk_key,
                        "chunk_size": len(chunk),
                        "status": status,
                    },
                )
                chunk_frames = []
                for identifier in chunk:
                    chunk_frames.extend(_fetch_single(identifier))
            else:
                raise
        except requests.RequestException as exc:
            if not _is_retryable_chunk_error(exc):
                raise
            logger.warning(
                "chunk_network_retry",
                extra={
                    "stage": "chunk_retry",
                    "chunk_key": chunk_key,
                    "chunk_size": len(chunk),
                    "error": exc.__class__.__name__,
                },
            )
            chunk_frames = []
            for identifier in chunk:
                chunk_frames.extend(_fetch_single(identifier))
        if chunk_frames:
            records.append(pd.concat(chunk_frames, ignore_index=True))
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
    fields: Sequence[str] | None = None,
    page_limit: int = 1000,
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
    effective_fields = tuple(fields) if fields else TESTITEM_QUERY_FIELDS
    max_limit = max(1, min(page_limit, 1000))
    effective_chunk = max(1, min(chunk_size, max_limit))
    base_params: list[tuple[str, str]] = [("format", "json"), ("limit", str(max_limit))]
    if effective_fields:
        base_params.append(("fields", ",".join(effective_fields)))
    query_string = urlencode(base_params)
    base = f"{cfg.chembl_base.rstrip('/')}/molecule.json?{query_string}"
    effective_timeout = timeout if timeout is not None else cfg.timeout_read
    for chunk in _chunked(valid, effective_chunk):
        for url_chunk in _split_chunk_for_url(chunk, base):
            records.extend(
                _fetch_testitem_chunk(
                    url_chunk,
                    base=base,
                    cfg=cfg,
                    client=client,
                    timeout=effective_timeout,
                )
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
    "TESTITEM_PUBCHEM_COLUMNS",
]
