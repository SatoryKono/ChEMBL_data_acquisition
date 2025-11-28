"""HTTP client for the Semantic Scholar API.

Changelog
========
* Extracted from :mod:`library.pubmed.query` to provide a dedicated client
  layer for Semantic Scholar interactions.
"""

from __future__ import annotations

import json
import logging
from typing import Any

import requests

from ..common.rate_limiter import RateLimiter
from ..config.models import RetryCfg, SemanticScholarCfg
from .base import RateLimitedRequestMixin

__all__ = [
    "fetch_semantic_scholar",
    "fetch_semantic_scholar_batch",
    "is_access_denied_error",
]

logger = logging.getLogger(__name__)

_REQUESTER = RateLimitedRequestMixin("semantic_scholar")

_SEMANTIC_SCHOLAR_FIELDS = "publicationTypes,externalIds,paperId,venue"
_SEMANTIC_SCHOLAR_HEADERS = {"Accept": "application/json"}

_ACCESS_DENIED_HINT = (
    "Semantic Scholar API access denied. Provide a valid API key via "
    "`sources.semantic_scholar.api_key` or disable the Semantic Scholar enrichment."
)


def is_access_denied_error(error: str | None) -> bool:
    """Return ``True`` when ``error`` represents an access denial response."""

    if not error:
        return False
    lowered = error.lower()
    return "http 401" in lowered or "http 403" in lowered or "forbidden" in lowered


def _format_error_message(error: str) -> str:
    """Attach actionable hints for well-known Semantic Scholar failures."""

    if is_access_denied_error(error):
        return f"{_ACCESS_DENIED_HINT} Original error: {error}"
    return error


def _build_headers(cfg: SemanticScholarCfg) -> dict[str, str]:
    """Return request headers for Semantic Scholar calls."""

    headers = dict(_SEMANTIC_SCHOLAR_HEADERS)
    if cfg.api_key:
        headers["x-api-key"] = cfg.api_key
    return headers


def fetch_semantic_scholar(
    session: requests.Session,
    pmid: str,
    cfg: SemanticScholarCfg | None = None,
    *,
    limiter: RateLimiter | None = None,
    retry_cfg: RetryCfg | None = None,
) -> dict[str, str]:
    """Retrieve Semantic Scholar metadata for ``pmid``."""

    cfg = cfg or SemanticScholarCfg()

    base = cfg.base.rstrip("/")
    url = f"{base}/paper/PMID:{pmid}"
    data, error = _REQUESTER.perform_request(
        session,
        url,
        cfg,
        limiter=limiter,
        retry_cfg=retry_cfg,
        headers=_build_headers(cfg),
        params={"fields": _SEMANTIC_SCHOLAR_FIELDS},
    )
    formatted_error = _format_error_message(error) if error else ""
    if error:
        if is_access_denied_error(error):
            logger.warning(
                "Semantic Scholar request failed for PMID %s due to access denial: %s",
                pmid,
                error,
            )
        else:
            logger.warning(
                "Semantic Scholar request failed for PMID %s with status: %s",
                pmid,
                error,
            )
    if error or not isinstance(data, dict):
        if not error:
            logger.warning(
                "Semantic Scholar request returned invalid response for PMID %s",
                pmid,
            )
        return {
            "scholar.PMID": pmid,
            "scholar.Venue": "",
            "scholar.PublicationTypes": "",
            "scholar.SemanticScholarId": "",
            "scholar.ExternalIds": "",
            "scholar.DOI": "",
            "scholar.Error": formatted_error or error or "Invalid response",
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
    cfg: SemanticScholarCfg | None = None,
    *,
    limiter: RateLimiter | None = None,
    retry_cfg: RetryCfg | None = None,
) -> list[dict[str, str]]:
    """Retrieve Semantic Scholar metadata for multiple PMIDs."""

    if not pmids:
        return []

    cfg = cfg or SemanticScholarCfg()

    base = cfg.base.rstrip("/")
    url = f"{base}/paper/batch"
    data, error = _REQUESTER.perform_request(
        session,
        url,
        cfg,
        limiter=limiter,
        retry_cfg=retry_cfg,
        headers=_build_headers(cfg),
        json={"ids": [f"PMID:{pmid}" for pmid in pmids]},
        method="POST",
        params={"fields": _SEMANTIC_SCHOLAR_FIELDS},
    )
    if error:
        if is_access_denied_error(error):
            logger.warning(
                "Semantic Scholar batch request failed for %d PMIDs due to access denial: %s",
                len(pmids),
                error,
            )
        else:
            for pmid in pmids:
                logger.warning(
                    "Semantic Scholar batch request failed for PMID %s with status: %s",
                    pmid,
                    error,
                )
        formatted_error = _format_error_message(error)
        return [
            {
                "scholar.PMID": pmid,
                "scholar.Venue": "",
                "scholar.PublicationTypes": "",
                "scholar.SemanticScholarId": "",
                "scholar.ExternalIds": "",
                "scholar.DOI": "",
                "scholar.Error": formatted_error,
            }
            for pmid in pmids
        ]

    if not isinstance(data, list):
        for pmid in pmids:
            logger.warning(
                "Semantic Scholar batch request returned invalid response for PMID %s",
                pmid,
            )
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
