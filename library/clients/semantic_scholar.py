"""HTTP client for the Semantic Scholar API.

Changelog
========
* Extracted from :mod:`library.pubmed.query` to provide a dedicated client
  layer for Semantic Scholar interactions.
"""

from __future__ import annotations

import json
from typing import Any

import requests

from ..config.models import RetryCfg, SemanticScholarCfg
from .pubmed import _do_request

__all__ = ["fetch_semantic_scholar", "fetch_semantic_scholar_batch"]

_SEMANTIC_SCHOLAR_FIELDS = "publicationTypes,externalIds,paperId,venue"
_SEMANTIC_SCHOLAR_HEADERS = {"Accept": "application/json"}


def fetch_semantic_scholar(
    session: requests.Session,
    pmid: str,
    sleep: float,
    cfg: SemanticScholarCfg | None = None,
    *,
    retry_cfg: RetryCfg | None = None,
) -> dict[str, str]:
    """Retrieve Semantic Scholar metadata for ``pmid``."""

    cfg = cfg or SemanticScholarCfg()
    base = cfg.base.rstrip("/")
    url = f"{base}/paper/PMID:{pmid}"
    timeout = (cfg.timeout_connect, cfg.timeout_read)
    data, error = _do_request(
        session,
        url,
        sleep * 5,
        headers=_SEMANTIC_SCHOLAR_HEADERS,
        params={"fields": _SEMANTIC_SCHOLAR_FIELDS},
        retries=cfg.retries,
        timeout=timeout,
        retry_cfg=retry_cfg,
    )
    if error or not isinstance(data, dict):
        return {
            "scholar.PMID": pmid,
            "scholar.Venue": "",
            "scholar.PublicationTypes": "",
            "scholar.SemanticScholarId": "",
            "scholar.ExternalIds": "",
            "scholar.DOI": "",
            "scholar.Error": error or "Invalid response",
        }

    external_ids = data.get("externalIds") or {}
    doi = external_ids.get("DOI") or ""
    pubtypes = data.get("publicationTypes") or []
    return {
        "scholar.PMID": pmid,
        "scholar.Venue": data.get("venue", ""),
        "scholar.PublicationTypes": "; ".join(pubtypes) if pubtypes else "",
        "scholar.SemanticScholarId": data.get("paperId", ""),
        "scholar.ExternalIds": json.dumps(external_ids, ensure_ascii=False),
        "scholar.DOI": doi,
        "scholar.Error": "",
    }


def fetch_semantic_scholar_batch(
    session: requests.Session,
    pmids: list[str],
    sleep: float,
    cfg: SemanticScholarCfg | None = None,
    *,
    retry_cfg: RetryCfg | None = None,
) -> list[dict[str, str]]:
    """Retrieve Semantic Scholar metadata for multiple PMIDs."""

    if not pmids:
        return []

    cfg = cfg or SemanticScholarCfg()
    base = cfg.base.rstrip("/")
    url = f"{base}/paper/batch"
    timeout = (cfg.timeout_connect, cfg.timeout_read)
    data, error = _do_request(
        session,
        url,
        sleep,
        headers=_SEMANTIC_SCHOLAR_HEADERS,
        json={"ids": [f"PMID:{pmid}" for pmid in pmids]},
        method="POST",
        params={"fields": _SEMANTIC_SCHOLAR_FIELDS},
        retries=cfg.retries,
        timeout=timeout,
        retry_cfg=retry_cfg,
    )
    if error:
        return [
            {
                "scholar.PMID": pmid,
                "scholar.Venue": "",
                "scholar.PublicationTypes": "",
                "scholar.SemanticScholarId": "",
                "scholar.ExternalIds": "",
                "scholar.DOI": "",
                "scholar.Error": error,
            }
            for pmid in pmids
        ]

    if not isinstance(data, list):
        return [
            {
                "scholar.PMID": pmid,
                "scholar.Venue": "",
                "scholar.PublicationTypes": "",
                "scholar.SemanticScholarId": "",
                "scholar.ExternalIds": "",
                "scholar.DOI": "",
                "scholar.Error": "Invalid batch response format",
            }
            for pmid in pmids
        ]

    def _resolve_pmid(item: dict[str, Any]) -> str | None:
        external_ids = item.get("externalIds")
        if isinstance(external_ids, dict):
            for key in ("PubMed", "PMID", "pubmed", "pubMed"):
                value = external_ids.get(key)
                if isinstance(value, str):
                    candidate = value.strip()
                    if candidate:
                        return candidate
                elif isinstance(value, list | tuple | set):
                    for entry in value:
                        if isinstance(entry, str):
                            candidate = entry.strip()
                            if candidate:
                                return candidate
                        elif entry is not None:
                            return str(entry)
                elif value is not None:
                    return str(value)
        for key in ("pmid", "PMID"):
            value = item.get(key)
            if isinstance(value, str):
                candidate = value.strip()
                if candidate:
                    return candidate
            elif value is not None:
                return str(value)
        return None

    lookup: dict[str, dict[str, Any]] = {}
    for entry in data:
        if not isinstance(entry, dict):
            continue
        pmid_value = _resolve_pmid(entry)
        if pmid_value:
            lookup.setdefault(pmid_value, entry)

    results: list[dict[str, str]] = []
    for pmid in pmids:
        item = lookup.get(pmid)
        if item is None:
            results.append(
                {
                    "scholar.PMID": pmid,
                    "scholar.Venue": "",
                    "scholar.PublicationTypes": "",
                    "scholar.SemanticScholarId": "",
                    "scholar.ExternalIds": "",
                    "scholar.DOI": "",
                    "scholar.Error": "Not found",
                }
            )
            continue

        external_ids = item.get("externalIds") or {}
        doi = external_ids.get("DOI") or ""
        pubtypes = item.get("publicationTypes") or []

        results.append(
            {
                "scholar.PMID": pmid,
                "scholar.Venue": item.get("venue", ""),
                "scholar.PublicationTypes": "; ".join(pubtypes) if pubtypes else "",
                "scholar.SemanticScholarId": item.get("paperId", ""),
                "scholar.ExternalIds": json.dumps(external_ids, ensure_ascii=False),
                "scholar.DOI": doi,
                "scholar.Error": "",
            }
        )

    return results
