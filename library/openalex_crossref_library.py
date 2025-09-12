"""OpenAlex and CrossRef query helpers.

These functions proxy to implementations in :mod:`library.pubmed_library` but
are exposed in a separate module to provide a clear separation of concerns.
"""

from __future__ import annotations

from urllib.parse import quote

import requests

from . import pubmed_library as _pl
from .config import CrossRefCfg, OpenAlexCfg
from .log import logger
from .rate_limiter import RateLimiter


def fetch_openalex(
    session: requests.Session,
    pmid: str,
    cfg: OpenAlexCfg,
    limiter: RateLimiter | None = None,
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
        Metadata returned by the OpenAlex API.

    Raises
    ------
    requests.RequestException
        Propagated if the HTTP request fails.

    """

    base = cfg.base.rstrip("/")
    url = f"{base}/works/pmid:{pmid}?mailto={quote(cfg.mailto)}"
    logger.info("request_start", extra={"stage": "request_start", "url": url})
    try:
        data = _pl.fetch_openalex(session, pmid, cfg=cfg, limiter=limiter)
    except requests.RequestException:
        logger.exception("request_fail", extra={"stage": "request_fail", "url": url})
        raise
    logger.info("request_ok", extra={"stage": "request_ok", "url": url})
    return data


def fetch_crossref(
    session: requests.Session,
    doi: str,
    cfg: CrossRefCfg,
    limiter: RateLimiter | None = None,
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
        Metadata returned by the CrossRef API.

    Raises
    ------
    requests.RequestException
        Propagated if the HTTP request fails.

    """

    base = cfg.base.rstrip("/")
    url = f"{base}/works/{quote(doi)}"
    logger.info("request_start", extra={"stage": "request_start", "url": url})
    try:
        data = _pl.fetch_crossref(session, doi, cfg=cfg, limiter=limiter)
    except requests.RequestException:
        logger.exception("request_fail", extra={"stage": "request_fail", "url": url})
        raise
    logger.info("request_ok", extra={"stage": "request_ok", "url": url})
    return data
