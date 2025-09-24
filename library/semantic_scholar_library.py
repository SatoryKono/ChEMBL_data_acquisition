"""Semantic Scholar query utilities.

This module wraps :func:`library.pubmed_library.fetch_semantic_scholar` to
provide a dedicated namespace for Semantic Scholar related functionality.
The underlying implementation already performs robust error handling and
request retries.
"""

from __future__ import annotations

import requests

from . import pubmed_library as _pl
from .config import SemanticScholarCfg


def fetch_semantic_scholar(
    session: requests.Session,
    pmid: str,
    sleep: float,
    cfg: SemanticScholarCfg | None = None,
) -> dict[str, str]:
    """Return Semantic Scholar metadata for ``pmid``.

    Parameters
    ----------
    session:
        :class:`requests.Session` instance used for the HTTP call.
    pmid:
        PubMed identifier of the article.
    sleep:
        Delay in seconds before making the request.  Semantic Scholar is more
        restrictive than PubMed; therefore the caller typically supplies a
        larger delay here.

    Returns
    -------
    dict
        Mapping of metadata fields to values.  Errors are encoded within the
        returned dictionary and never raise exceptions.

    """
    return _pl.fetch_semantic_scholar(session, pmid, sleep, cfg=cfg)


def fetch_semantic_scholar_batch(
    session: requests.Session,
    pmids: list[str],
    sleep: float,
    cfg: SemanticScholarCfg | None = None,
) -> list[dict[str, str]]:
    """Return Semantic Scholar metadata for a batch of ``pmids``.

    Parameters
    ----------
    session:
        :class:`requests.Session` instance used for the HTTP call.
    pmids:
        A list of PubMed identifiers.
    sleep:
        Delay in seconds before making the request.

    Returns
    -------
    list of dict
        A list of metadata mappings. Errors are encoded within each
        dictionary and never raise exceptions.

    """
    return _pl.fetch_semantic_scholar_batch(session, pmids, sleep, cfg=cfg)
