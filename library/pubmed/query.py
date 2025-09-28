"""HTTP queries for publication metadata.

This module provides functions for retrieving metadata from PubMed,
Semantic Scholar, OpenAlex and Crossref.  It also contains helpers for
issuing robust HTTP requests.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any
from urllib.parse import quote
from xml.etree import ElementTree as ET

import pandas as pd
import requests

from ..config import CrossRefCfg, OpenAlexCfg, PubMedCfg, SemanticScholarCfg
from ..utils.logging import logger
from ..rate_limiter import get_limiter, sleep
from .parsing import EMPTY_PUBMED, combine, parse_pubmed_article

if TYPE_CHECKING:  # pragma: no cover - imported for type checking only
    from ..rate_limiter import RateLimiter


__all__ = [
    "read_pmids",
    "_make_request",
    "_handle_response",
    "_do_request",
    "fetch_pubmed_batch",
    "fetch_pubmed",
    "fetch_semantic_scholar",
    "fetch_semantic_scholar_batch",
    "fetch_openalex",
    "fetch_crossref",
]


_SEMANTIC_SCHOLAR_FIELDS = "publicationTypes,externalIds,paperId,venue"
_SEMANTIC_SCHOLAR_HEADERS = {"Accept": "application/json"}
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


def _make_request(
    session: requests.Session,
    url: str,
    expect_json: bool,
    method: str,
    timeout: float | tuple[float, float],
    **kwargs: Any,
) -> tuple[int, str, dict[str, Any] | str | None, str]:
    """Issue a single HTTP request and parse its body."""

    if method.upper() == "POST":
        request: Callable[..., requests.Response] = session.post
    else:
        request = session.get
    with request(url, timeout=timeout, **kwargs) as resp:
        status_code = resp.status_code
        text = resp.text
        if expect_json:
            try:
                content = resp.json()
            except ValueError as exc:
                return status_code, text, None, str(exc)
            return status_code, text, content, ""
        return status_code, text, text or "", ""


def _handle_response(
    url: str,
    status_code: int,
    text: str,
    content: dict[str, Any] | str | None,
    parse_error: str,
    expect_json: bool,
    attempt: int,
    retries: int,
) -> tuple[dict[str, Any] | str | None, str, bool]:
    """Process an HTTP response and decide whether to retry."""

    if status_code in (429, 500, 502, 503, 504):
        if attempt >= retries:
            logger.info(
                "request_fail",
                extra={"stage": "request_fail", "url": url, "status": status_code},
            )
            return None, f"HTTP {status_code}: {text[:100]}", False
        return None, "", True
    if status_code == 404:
        logger.info(
            "request_fail",
            extra={"stage": "request_fail", "url": url, "status": status_code},
        )
        return None, "PMID not found", False
    if status_code == 400:
        logger.info(
            "request_fail",
            extra={"stage": "request_fail", "url": url, "status": status_code},
        )
        return None, f"Bad request: {text[:100]}", False
    if status_code != 200:
        logger.info(
            "request_fail",
            extra={"stage": "request_fail", "url": url, "status": status_code},
        )
        return None, f"HTTP {status_code}: {text[:100]}", False
    if expect_json:
        if parse_error:
            logger.info("request_fail", extra={"stage": "request_fail", "url": url})
            return None, f"Invalid JSON: {parse_error}", False
        logger.info(
            "request_ok",
            extra={"stage": "request_ok", "url": url, "status": status_code},
        )
        return content, "", False
    logger.info(
        "request_ok",
        extra={"stage": "request_ok", "url": url, "status": status_code},
    )
    return content or "", "", False


def _do_request(
    session: requests.Session,
    url: str,
    delay: float,
    expect_json: bool = True,
    retries: int = 2,
    method: str = "GET",
    timeout: float | tuple[float, float] = 10,
    **kwargs: Any,
) -> tuple[dict[str, Any] | str | None, str]:
    """Perform an HTTP request with retry logic."""

    for attempt in range(retries + 1):
        event = "request_start" if attempt == 0 else "request_retry"
        logger.info(event, extra={"stage": event, "url": url, "attempt": attempt + 1})
        if attempt:
            sleep(delay * attempt)

        try:
            status_code, text, content, parse_error = _make_request(
                session, url, expect_json, method, timeout, **kwargs
            )
        except requests.RequestException as exc:
            if attempt >= retries:  # pragma: no cover - network errors
                logger.exception(
                    "request_fail", extra={"stage": "request_fail", "url": url}
                )
                return None, str(exc)
            continue

        data, error, retry = _handle_response(
            url, status_code, text, content, parse_error, expect_json, attempt, retries
        )
        if retry:
            continue
        return data, error

    logger.info("request_fail", extra={"stage": "request_fail", "url": url})
    return None, "Request failed"


def fetch_pubmed_batch(
    session: requests.Session,
    pmids: list[str],
    sleep: float,
    cfg: PubMedCfg | None = None,
) -> list[dict[str, str]]:
    """Fetch metadata for multiple PMIDs using a single API request."""

    cfg = cfg or PubMedCfg()
    ids = ",".join(pmids)
    base = cfg.base.rstrip("/")
    url = f"{base}/efetch.fcgi?db=pubmed&id={ids}&retmode=xml"
    timeout = (cfg.timeout_connect, cfg.timeout_read)
    text, error = _do_request(
        session, url, sleep, expect_json=False, retries=cfg.retries, timeout=timeout
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
        root = ET.fromstring(text)  # type: ignore[arg-type]
    except ET.ParseError as exc:  # pragma: no cover - malformed XML
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
    session: requests.Session, pmid: str, sleep: float, cfg: PubMedCfg | None = None
) -> dict[str, str]:
    """Fetch metadata for a single PMID."""

    return fetch_pubmed_batch(session, [pmid], sleep, cfg=cfg)[0]


def fetch_semantic_scholar(
    session: requests.Session,
    pmid: str,
    sleep: float,
    cfg: SemanticScholarCfg | None = None,
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


def fetch_openalex(
    session: requests.Session,
    pmid: str,
    *,
    cfg: OpenAlexCfg,
    limiter: RateLimiter | None = None,
) -> dict[str, str]:
    """Retrieve OpenAlex metadata for ``pmid``."""

    if limiter is None:
        limiter = get_limiter("openalex", cfg.rps, cfg.burst)
    limiter.acquire()
    delay = 1 / cfg.rps if cfg.rps else 0
    base = cfg.base.rstrip("/")
    url = f"{base}/works/pmid:{pmid}?mailto={quote(cfg.mailto)}"
    timeout = (cfg.timeout_connect, cfg.timeout_read)
    data, error = _do_request(session, url, delay, retries=cfg.retries, timeout=timeout)

    if error or not isinstance(data, dict):
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
    mesh_entries = data.get("mesh") or []
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
        "OpenAlex.PublicationTypes": data.get("type", ""),
        "OpenAlex.TypeCrossref": data.get("type_crossref", ""),
        "OpenAlex.Genre": data.get("genre", ""),
        "OpenAlex.Id": data.get("id", ""),
        "OpenAlex.Venue": data.get("host_venue", {}).get("display_name", ""),
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

    if limiter is None:
        limiter = get_limiter("crossref", cfg.rps, cfg.burst)
    limiter.acquire()
    delay = 1 / cfg.rps if cfg.rps else 0
    base = cfg.base.rstrip("/")
    url = f"{base}/works/{quote(doi, safe='')}?mailto={quote(cfg.mailto)}"
    timeout = (cfg.timeout_connect, cfg.timeout_read)
    data, error = _do_request(session, url, delay, retries=cfg.retries, timeout=timeout)

    if error or not isinstance(data, dict):
        return {
            "crossref.Type": "",
            "crossref.Subtype": "",
            "crossref.Title": "",
            "crossref.Subtitle": "",
            "crossref.Subject": "",
            "crossref.Error": error or "Invalid response",
        }
    message = data.get("message", {})
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
