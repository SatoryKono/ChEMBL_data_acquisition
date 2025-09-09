"""OpenAlex and CrossRef query helpers.

These functions proxy to implementations in :mod:`library.pubmed_library` but
are exposed in a separate module to provide a clear separation of concerns.
"""

from __future__ import annotations

from typing import Dict
import logging

import requests

from . import pubmed_library as _pl

logger = logging.getLogger(__name__)


def fetch_openalex(
    session: requests.Session,
    pmid: str,
    sleep: float,
    *,
    timeout: float = _pl.TIMEOUT,
) -> Dict[str, str]:
    """Return OpenAlex metadata for ``pmid``.

    Parameters
    ----------
    session: requests.Session
        Session used for the HTTP request.
    pmid: str
        PubMed identifier.
    sleep: float
        Delay before making the request in seconds.
    timeout: float, optional
        Maximum seconds to wait for the HTTP response.

    """
    return _pl.fetch_openalex(session, pmid, sleep, timeout=timeout)


def fetch_crossref(
    session: requests.Session,
    doi: str,
    sleep: float,
    *,
    timeout: float = _pl.TIMEOUT,
) -> Dict[str, str]:
    """Return CrossRef metadata for ``doi``.

    Parameters
    ----------
    session: requests.Session
        Session used for the HTTP request.
    doi: str
        Digital Object Identifier of the article.
    sleep: float
        Delay before making the request in seconds.
    timeout: float, optional
        Maximum seconds to wait for the HTTP response.

    """
    return _pl.fetch_crossref(session, doi, sleep, timeout=timeout)
