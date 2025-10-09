"""OpenAlex and CrossRef query helpers.

These functions proxy to implementations in :mod:`library.integration.pubmed_library` but
are exposed in a separate module to provide a clear separation of concerns.
"""

from __future__ import annotations

from typing import Any, Callable, Iterable, cast

import requests

from ..clients import crossref as crossref_client
from ..clients import openalex as openalex_client
from ..config import CrossRefCfg, OpenAlexCfg, RetryCfg
from ..pubmed import query as pubmed_query
from ..common.rate_limiter import RateLimiter


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
    if error or not isinstance(raw, dict):
        return {
            "OpenAlex.PublicationTypes": "",
            "OpenAlex.TypeCrossref": "",
            "OpenAlex.Genre": "",
            "OpenAlex.Id": "",
            "OpenAlex.Venue": "",
            "OpenAlex.MeshDescriptors": "",
            "OpenAlex.MeshQualifiers": "",
            "OpenAlex.Error": error or "Invalid response",
        }

    mesh_entries = raw.get("mesh") or []
    descriptors: list[str] = []
    qualifiers: list[str] = []
    for entry in mesh_entries:
        descriptor = entry.get("descriptor_name")
        if descriptor:
            descriptors.append(descriptor)
        for qualifier in entry.get("qualifiers") or []:
            qualifier_name = qualifier.get("qualifier_name")
            if qualifier_name:
                qualifiers.append(qualifier_name)

    return {
        "OpenAlex.PublicationTypes": raw.get("type", ""),
        "OpenAlex.TypeCrossref": raw.get("type_crossref", ""),
        "OpenAlex.Genre": raw.get("genre", ""),
        "OpenAlex.Id": raw.get("id", ""),
        "OpenAlex.Venue": raw.get("host_venue", {}).get("display_name", ""),
        "OpenAlex.MeshDescriptors": _combine_mesh(descriptors),
        "OpenAlex.MeshQualifiers": _combine_mesh(qualifiers),
        "OpenAlex.Error": "",
    }


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
        return {
            "crossref.Type": "",
            "crossref.Subtype": "",
            "crossref.Title": "",
            "crossref.Subtitle": "",
            "crossref.Subject": "",
            "crossref.Error": "Missing DOI",
        }

    raw, error = crossref_client.fetch_crossref(
        session,
        doi,
        cfg=cfg,
        limiter=limiter,
        retry_cfg=retry_cfg,
    )
    if error or not isinstance(raw, dict):
        return {
            "crossref.Type": "",
            "crossref.Subtype": "",
            "crossref.Title": "",
            "crossref.Subtitle": "",
            "crossref.Subject": "",
            "crossref.Error": error or "Invalid response",
        }

    message: dict[str, Any] = raw.get("message", {})
    title = message.get("title") or [""]
    subtitle = message.get("subtitle") or [""]
    subject = "; ".join(message.get("subject") or [])

    return {
        "crossref.Type": message.get("type", ""),
        "crossref.Subtype": message.get("subtype", ""),
        "crossref.Title": title[0] if title else "",
        "crossref.Subtitle": subtitle[0] if subtitle else "",
        "crossref.Subject": subject,
        "crossref.Error": "",
    }


_COMBINE: Callable[[Iterable[str]], str] = cast(
    Callable[[Iterable[str]], str], getattr(pubmed_query, "combine")
)


def _combine_mesh(tokens: Iterable[str]) -> str:
    """Proxy to :func:`library.pubmed.query.combine` with precise typing."""

    return _COMBINE(tokens)
