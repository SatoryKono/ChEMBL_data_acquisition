"""Semantic Scholar query utilities.

This module provides a thin processing wrapper around
``library.clients.semantic_scholar``.  The client encapsulates HTTP concerns,
while this module offers a semantic namespace for higher-level pipelines.
"""

from __future__ import annotations

import requests

from ..clients import semantic_scholar as _client
from ..config import RetryCfg, SemanticScholarCfg


def fetch_semantic_scholar(
    session: requests.Session,
    pmid: str,
    sleep: float,
    cfg: SemanticScholarCfg | None = None,
    *,
    retry_cfg: RetryCfg | None = None,
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
    return _client.fetch_semantic_scholar(
        session,
        pmid,
        sleep,
        cfg=cfg,
        retry_cfg=retry_cfg,
    )


def fetch_semantic_scholar_batch(
    session: requests.Session,
    pmids: list[str],
    sleep: float,
    cfg: SemanticScholarCfg | None = None,
    *,
    retry_cfg: RetryCfg | None = None,
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
    return _client.fetch_semantic_scholar_batch(
        session,
        pmids,
        sleep,
        cfg=cfg,
        retry_cfg=retry_cfg,
    )
