"""Document helpers for the ChEMBL API."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from time import monotonic
from typing import Any

import pandas as pd
import requests

from library.clients import ChemblClient, _chunked

from ...common.log import logger
from ...config import ApiCfg

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

INVALID_DOCUMENT_IDS = {"", "#N/A"}

MAX_DOCUMENT_CHUNK_SIZE = 20


def get_documents(
    ids: Iterable[str],
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    chunk_size: int = 5,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Fetch document records for ``ids`` from the ChEMBL API.

    Parameters
    ----------
    ids:
        ChEMBL document identifiers to retrieve.
    cfg:
        API configuration providing base URL and timeouts.
    client:
        :class:`ChemblClient` instance used for HTTP requests and caching.
    chunk_size:
        Maximum number of identifiers per HTTP request.
    timeout:
        Optional override for the read timeout in seconds.

    Returns
    -------
    pandas.DataFrame
        Data frame containing the retrieved document metadata. Missing
        identifiers result in an empty frame.
    """
    base = f"{cfg.chembl_base.rstrip('/')}/document.json?format=json"
    effective_timeout = timeout if timeout is not None else cfg.timeout_read
    seen: set[str] = set()
    unique_ids: list[str] = []

    for identifier in ids:
        if identifier in INVALID_DOCUMENT_IDS:
            continue
        if identifier in seen:
            continue
        seen.add(identifier)
        unique_ids.append(identifier)

    if not unique_ids:
        return pd.DataFrame(columns=DOCUMENT_COLUMNS)

    records: list[dict[str, Any]] = []
    effective_chunk_size = max(1, min(chunk_size, MAX_DOCUMENT_CHUNK_SIZE))
    if effective_chunk_size != chunk_size:
        logger.info(
            "chembl_document_chunk_size_capped",
            extra={
                "requested": chunk_size,
                "applied": effective_chunk_size,
                "max": MAX_DOCUMENT_CHUNK_SIZE,
            },
        )

    for chunk in _chunked(unique_ids, effective_chunk_size):
        records.extend(
            _fetch_documents_chunk(
                chunk,
                base_url=base,
                cfg=cfg,
                client=client,
                timeout=effective_timeout,
            )
        )

    if not records:
        return pd.DataFrame(columns=DOCUMENT_COLUMNS)

    df = pd.DataFrame(records)
    return df.reindex(columns=DOCUMENT_COLUMNS)


def _fetch_documents_chunk(
    chunk_ids: Sequence[str],
    *,
    base_url: str,
    cfg: ApiCfg,
    client: ChemblClient,
    timeout: float,
) -> list[dict[str, Any]]:
    """Return document records for ``chunk_ids`` with timeout-aware fallback."""

    if not chunk_ids:
        return []

    url = f"{base_url}&document_chembl_id__in={','.join(chunk_ids)}"
    start_time = monotonic()
    try:
        data = client.request_json(url, cfg=cfg, timeout=timeout)
    except requests.ReadTimeout as exc:
        if len(chunk_ids) <= 1:
            raise
        elapsed = monotonic() - start_time
        logger.warning(
            "chembl_document_timeout_split",
            extra={
                "chunk_size": len(chunk_ids),
                "ids": list(chunk_ids),
                "timeout": timeout,
                "error": str(exc),
                "elapsed": elapsed,
            },
        )
        midpoint = max(1, len(chunk_ids) // 2)
        left = chunk_ids[:midpoint]
        right = chunk_ids[midpoint:]
        records: list[dict[str, Any]] = []
        if left:
            records.extend(
                _fetch_documents_chunk(
                    left,
                    base_url=base_url,
                    cfg=cfg,
                    client=client,
                    timeout=timeout,
                )
            )
        if right:
            records.extend(
                _fetch_documents_chunk(
                    right,
                    base_url=base_url,
                    cfg=cfg,
                    client=client,
                    timeout=timeout,
                )
            )
        return records
    except requests.RequestException:
        raise

    items = data.get("documents") or data.get("document") or []
    records = [
        _normalise_document_record(item) for item in items if isinstance(item, Mapping)
    ]
    return [record for record in records if record]


def _normalise_document_record(item: Mapping[str, Any]) -> dict[str, Any]:
    """Return a normalised document record extracted from ``item``."""

    return {
        "document_chembl_id": item.get("document_chembl_id"),
        "title": item.get("title"),
        "abstract": item.get("abstract"),
        "doi": item.get("doi"),
        "year": item.get("year"),
        "journal": item.get("journal_full_title"),
        "journal_abbrev": item.get("journal"),
        "volume": item.get("volume"),
        "issue": item.get("issue"),
        "first_page": item.get("first_page"),
        "last_page": item.get("last_page"),
        "pubmed_id": item.get("pubmed_id"),
        "authors": item.get("authors"),
        "source": "ChEMBL",
    }


__all__ = ["get_documents"]
