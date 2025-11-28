"""OpenAlex and CrossRef query helpers.

These functions proxy to implementations in :mod:`library.pubmed.integration.pubmed_library` but
are exposed in a separate module to provide a clear separation of concerns.
"""

from __future__ import annotations

import requests

from library.crossref.clients import crossref as crossref_client
from ..clients import openalex as openalex_client
from ..common.rate_limiter import RateLimiter
from ..config import CrossRefCfg, OpenAlexCfg, RetryCfg
from ._normalizers import (
    normalize_crossref_response,
    normalize_openalex_response,
)


def fetch_openalex(
    session: requests.Session,
    pmid: str,
    cfg: OpenAlexCfg,
    limiter: RateLimiter | None = None,
    *,
    retry_cfg: RetryCfg | None = None,
) -> dict[str, str]:
    """Return OpenAlex metadata for ``pmid``.

    Parameters
    ----------
    session: requests.Session
        Session used for the HTTP request.
    pmid: str
        PubMed identifier.
    cfg: OpenAlexCfg
        Configuration specifying base URL, timeouts and rate limits.
    limiter: RateLimiter | None
        Shared rate limiter controlling request throughput.


    Returns
    -------
    dict
        Normalized metadata returned by the OpenAlex API.

    """

    raw, error = openalex_client.fetch_openalex(
        session,
        pmid,
        cfg=cfg,
        limiter=limiter,
        retry_cfg=retry_cfg,
    )
    return normalize_openalex_response(raw, error)


def fetch_crossref(
    session: requests.Session,
    doi: str,
    cfg: CrossRefCfg,
    limiter: RateLimiter | None = None,
    *,
    retry_cfg: RetryCfg | None = None,
) -> dict[str, str]:
    """Return CrossRef metadata for ``doi``.

    Parameters
    ----------
    session: requests.Session
        Session used for the HTTP request.
    doi: str
        Digital Object Identifier of the article.
    cfg: CrossRefCfg
        Configuration specifying base URL, timeouts and rate limits.
    limiter: RateLimiter | None
        Shared rate limiter controlling request throughput.


    Returns
    -------
    dict
        Normalized metadata returned by the CrossRef API.

    """

    if not doi:
        return normalize_crossref_response(None, "Missing DOI")

    raw, error = crossref_client.fetch_crossref(
        session,
        doi,
        cfg=cfg,
        limiter=limiter,
        retry_cfg=retry_cfg,
    )
    return normalize_crossref_response(raw, error)
