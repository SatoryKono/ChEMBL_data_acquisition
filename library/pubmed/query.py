"""HTTP queries for publication metadata.

This module provides functions for retrieving metadata from PubMed,
OpenAlex and Crossref.  It also contains helpers for
issuing robust HTTP requests.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING
from xml.etree import ElementTree as ET

import pandas as pd
import requests

from ..clients import (
    crossref as crossref_client,
)
from ..clients import (
    openalex as openalex_client,
)
from ..clients import (
    pubmed as pubmed_client,
)
from ..clients import (
    semantic_scholar as semantic_client,
)
from ..config import (
    CrossRefCfg,
    OpenAlexCfg,
    PubMedCfg,
    RetryCfg,
    SemanticScholarCfg,
)
from .parsing import EMPTY_PUBMED, combine, parse_pubmed_article

if TYPE_CHECKING:  # pragma: no cover - imported for type checking only
    from ..common.rate_limiter import RateLimiter


__all__ = [
    "read_pmids",
    "fetch_pubmed_batch",
    "fetch_pubmed",
    "fetch_semantic_scholar",
    "fetch_semantic_scholar_batch",
    "fetch_openalex",
    "fetch_crossref",
]


_PUBMED_FALLBACK_ENCODINGS: tuple[str, ...] = (
    "utf-8-sig",
    "cp1251",
    "latin1",
)


def read_pmids(path: str | Path, cfg: PubMedCfg | None = None) -> pd.DataFrame:
    """Read the ``PMID`` column from ``path`` as a DataFrame.

    Parameters
    ----------
    path:
        CSV file containing a ``PMID`` column.
    cfg:
        Optional :class:`PubMedCfg` providing fallback encodings.

    Returns
    -------
    pandas.DataFrame
        DataFrame with a single ``PMID`` column containing non-empty values.

    Raises
    ------
    ValueError
        If the file cannot be decoded using the configured encodings or
        the ``PMID`` column is missing.
    """

    path = Path(path)
    last_exc: Exception | None = None
    encodings = getattr(cfg or PubMedCfg(), "encodings", _PUBMED_FALLBACK_ENCODINGS)
    for enc in encodings:
        try:
            df = pd.read_csv(path, encoding=enc, dtype=str)
        except UnicodeDecodeError as exc:  # pragma: no cover - depends on filesystem
            last_exc = exc
            continue
        if "PMID" not in df.columns:
            raise ValueError("Input CSV must contain 'PMID' column")
        pmid_df = df[["PMID"]].copy()
        pmid_df["PMID"] = pmid_df["PMID"].fillna("").astype(str).str.strip()
        pmid_df = pmid_df[pmid_df["PMID"].astype(bool)]
        return pmid_df
    raise ValueError(
        f"Could not decode {path} with encodings {encodings}. Last error: {last_exc}"
    )


def fetch_pubmed_batch(
    session: requests.Session,
    pmids: list[str],
    sleep: float,
    cfg: PubMedCfg | None = None,
    *,
    client: pubmed_client.PubMedClient | None = None,
    retry_cfg: RetryCfg | None = None,
) -> list[dict[str, str]]:
    """Fetch and parse PubMed metadata for multiple PMIDs."""

    if client is not None and cfg is not None:
        raise ValueError("Specify either cfg or client, not both")

    resolved_client = client or pubmed_client.PubMedClient(cfg=cfg)
    xml_text, error = resolved_client.fetch_pubmed_batch(
        session, pmids, sleep, retry_cfg=retry_cfg
    )
    results: list[dict[str, str]] = []
    if error:
        for pid in pmids:
            res = EMPTY_PUBMED.copy()
            res["PubMed.PMID"] = pid
            res["PubMed.Error"] = error
            results.append(res)
        return results
    try:
        if xml_text is None:
            raise ValueError("Empty PubMed response")
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:  # pragma: no cover - malformed XML
        for pid in pmids:
            res = EMPTY_PUBMED.copy()
            res["PubMed.PMID"] = pid
            res["PubMed.Error"] = str(exc)
            results.append(res)
        return results
    except ValueError as exc:
        for pid in pmids:
            res = EMPTY_PUBMED.copy()
            res["PubMed.PMID"] = pid
            res["PubMed.Error"] = str(exc)
            results.append(res)
        return results
    articles = root.findall(".//PubmedArticle")
    parsed: dict[str, dict[str, str]] = {}
    for art in articles:
        rec = EMPTY_PUBMED.copy()
        rec.update(parse_pubmed_article(art))
        parsed[rec.get("PubMed.PMID", "")] = rec
    for pid in pmids:
        if pid in parsed:
            results.append(parsed[pid])
        else:
            missing = EMPTY_PUBMED.copy()
            missing["PubMed.PMID"] = pid
            missing["PubMed.Error"] = "No PubmedArticle"
            results.append(missing)
    return results


def fetch_pubmed(
    session: requests.Session,
    pmid: str,
    sleep: float,
    cfg: PubMedCfg | None = None,
    *,
    client: pubmed_client.PubMedClient | None = None,
    retry_cfg: RetryCfg | None = None,
) -> dict[str, str]:
    """Fetch metadata for a single PMID."""

    resolved_client = client or (pubmed_client.PubMedClient(cfg=cfg))
    return fetch_pubmed_batch(
        session,
        [pmid],
        sleep,
        client=resolved_client,
        retry_cfg=retry_cfg,
    )[0]


def fetch_semantic_scholar(
    session: requests.Session,
    pmid: str,
    sleep: float,
    cfg: SemanticScholarCfg | None = None,
    *,
    limiter: "RateLimiter" | None = None,
    retry_cfg: RetryCfg | None = None,
) -> dict[str, str]:
    """Retrieve Semantic Scholar metadata for ``pmid``."""
    return semantic_client.fetch_semantic_scholar(
        session,
        pmid,
        sleep,
        cfg=cfg,
        limiter=limiter,
        retry_cfg=retry_cfg,
    )


def fetch_semantic_scholar_batch(
    session: requests.Session,
    pmids: list[str],
    sleep: float,
    cfg: SemanticScholarCfg | None = None,
    *,
    limiter: "RateLimiter" | None = None,
    retry_cfg: RetryCfg | None = None,
) -> list[dict[str, str]]:
    """Retrieve Semantic Scholar metadata for multiple PMIDs."""
    return semantic_client.fetch_semantic_scholar_batch(
        session,
        pmids,
        sleep,
        cfg=cfg,
        limiter=limiter,
        retry_cfg=retry_cfg,
    )


def fetch_openalex(
    session: requests.Session,
    pmid: str,
    *,
    cfg: OpenAlexCfg,
    limiter: RateLimiter | None = None,
    retry_cfg: RetryCfg | None = None,
) -> dict[str, str]:
    """Retrieve OpenAlex metadata for ``pmid``."""
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
        d = entry.get("descriptor_name")
        if d:
            descriptors.append(d)
        for q in entry.get("qualifiers") or []:
            qn = q.get("qualifier_name")
            if qn:
                qualifiers.append(qn)
    return {
        "OpenAlex.PublicationTypes": raw.get("type", ""),
        "OpenAlex.TypeCrossref": raw.get("type_crossref", ""),
        "OpenAlex.Genre": raw.get("genre", ""),
        "OpenAlex.Id": raw.get("id", ""),
        "OpenAlex.Venue": raw.get("host_venue", {}).get("display_name", ""),
        "OpenAlex.MeshDescriptors": combine(descriptors),
        "OpenAlex.MeshQualifiers": combine(qualifiers),
        "OpenAlex.Error": "",
    }


def fetch_crossref(
    session: requests.Session,
    doi: str,
    *,
    cfg: CrossRefCfg,
    limiter: RateLimiter | None = None,
    retry_cfg: RetryCfg | None = None,
) -> dict[str, str]:
    """Retrieve Crossref metadata for a given DOI."""

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
    message = raw.get("message", {})
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
