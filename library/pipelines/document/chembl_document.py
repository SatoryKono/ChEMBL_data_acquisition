"""Document helpers for the ChEMBL API."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import pandas as pd

from library.clients import ChemblClient
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
    records: list[dict[str, Any]] = []
    seen: set[str] = set()
    chunk: list[str] = []

    def fetch_chunk(chunk_ids: list[str]) -> None:
        url = f"{base}&document_chembl_id__in={','.join(chunk_ids)}"
        data = client.request_json(url, cfg=cfg, timeout=effective_timeout)
        items = data.get("documents") or data.get("document") or []
        for item in items:
            record = {
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
            records.append(record)

    for identifier in ids:
        if identifier in INVALID_DOCUMENT_IDS:
            continue
        if identifier in seen:
            continue
        seen.add(identifier)
        chunk.append(identifier)
        if len(chunk) >= chunk_size:
            fetch_chunk(chunk)
            chunk = []

    if chunk:
        fetch_chunk(chunk)

    if not seen or not records:
        return pd.DataFrame(columns=DOCUMENT_COLUMNS)

    df = pd.DataFrame(records)
    return df.reindex(columns=DOCUMENT_COLUMNS)


__all__ = ["get_documents"]
