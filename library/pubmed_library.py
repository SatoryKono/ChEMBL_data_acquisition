"""Utilities for retrieving and merging publication metadata.

The functions in this module fetch information from PubMed, Semantic Scholar,
OpenAlex and Crossref and consolidate it into tabular form.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from collections.abc import Callable, Iterable, Sequence
from datetime import date
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
)
from urllib.parse import quote
from xml.etree import ElementTree as ET

import pandas as pd
import requests

from .config import (
    Config,
    CrossRefCfg,
    OpenAlexCfg,
    PubMedCfg,
    SemanticScholarCfg,
    session_with_retry,
)
from .csv_utils import write_csv_deterministic
from .log import logger
from .rate_limiter import get_limiter, sleep

if TYPE_CHECKING:
    from .rate_limiter import RateLimiter


def read_pmids(path: str | Path, cfg: PubMedCfg | None = None) -> list[str]:
    """Read PMID column from a CSV file.

    Parameters
    ----------
    path:
        CSV file containing a ``PMID`` column.
    cfg:
        Optional :class:`PubMedCfg` providing fallback encodings.

    Returns
    -------
    list of str
        Extracted PMIDs.

    """
    path = Path(path)
    last_exc: Exception | None = None
    encodings = (cfg or PubMedCfg()).encodings
    for enc in encodings:
        try:
            with path.open(encoding=enc, newline="") as f:
                reader = csv.DictReader(f)
                if reader.fieldnames is None or "PMID" not in reader.fieldnames:
                    raise ValueError("Input CSV must contain 'PMID' column")
                return [pmid for row in reader if (pmid := row.get("PMID", "").strip())]
        except UnicodeDecodeError as exc:
            last_exc = exc
            continue
    raise ValueError(
        f"Could not decode {path} with encodings {encodings}. Last error: {last_exc}"
    )


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
    """Perform an HTTP request with retry and error handling.

    Parameters
    ----------
    session:
        Requests session used to perform the call.
    url:
        Endpoint to query.
    delay:
        Base number of seconds to wait between attempts.
    expect_json:
        Whether to parse the response as JSON.
    retries:
        Number of additional attempts after the initial one.
    method:
        HTTP method to use, either "GET" or "POST".
    timeout:

        Maximum seconds to wait for each HTTP request. May be a single float
        or a ``(connect, read)`` tuple.

    **kwargs:
        Additional parameters passed to ``session.get`` or ``session.post``.

    Returns
    -------
    tuple
        Pair of parsed data (dict/str/``None``) and an error message. The error
        string is empty on success.

    """
    for attempt in range(retries + 1):
        event = "request_start" if attempt == 0 else "request_retry"
        logger.info(event, extra={"stage": event, "url": url, "attempt": attempt + 1})
        if attempt:
            sleep(delay * attempt)

        try:
            if method.upper() == "POST":
                request: Callable[..., requests.Response] = session.post
            else:
                request = session.get
            with request(url, timeout=timeout, **kwargs) as resp:
                status_code = resp.status_code
                text = resp.text
                try:
                    content = resp.json() if expect_json else text
                except ValueError as exc:
                    content = None
                    parse_error = str(exc)
                else:
                    parse_error = ""
        except requests.RequestException as exc:
            if attempt >= retries:  # pragma: no cover - network errors
                logger.exception(
                    "request_fail", extra={"stage": "request_fail", "url": url}
                )
                return None, str(exc)
            continue
        if status_code in (429, 500, 502, 503, 504):
            if attempt >= retries:
                logger.info(
                    "request_fail",
                    extra={
                        "stage": "request_fail",
                        "url": url,
                        "status": status_code,
                    },
                )
                return None, f"HTTP {status_code}: {text[:100]}"
            continue
        if status_code == 404:
            logger.info(
                "request_fail",
                extra={"stage": "request_fail", "url": url, "status": status_code},
            )
            return None, "PMID not found"
        if status_code == 400:
            logger.info(
                "request_fail",
                extra={"stage": "request_fail", "url": url, "status": status_code},
            )
            return None, f"Bad request: {text[:100]}"
        if status_code != 200:
            logger.info(
                "request_fail",
                extra={"stage": "request_fail", "url": url, "status": status_code},
            )
            return None, f"HTTP {status_code}: {text[:100]}"
        if expect_json:
            if parse_error:
                logger.info("request_fail", extra={"stage": "request_fail", "url": url})
                return None, f"Invalid JSON: {parse_error}"
            logger.info(
                "request_ok",
                extra={"stage": "request_ok", "url": url, "status": status_code},
            )
            return content, ""
        logger.info(
            "request_ok",
            extra={"stage": "request_ok", "url": url, "status": status_code},
        )
        return content or "", ""
    logger.info("request_fail", extra={"stage": "request_fail", "url": url})
    return None, "Request failed"


def fetch_pubmed_batch(
    session: requests.Session,
    pmids: list[str],
    sleep: float,
    cfg: PubMedCfg | None = None,
) -> list[dict[str, str]]:
    """Fetch metadata for multiple PMIDs using a single API request.

    Parameters
    ----------
    session:
        Active :class:`requests.Session`.
    pmids:
        List of PubMed identifiers.
    sleep:
        Seconds to pause after the request.
    cfg:
        Optional :class:`PubMedCfg` with API settings.

    Returns
    -------
    list of dict
        One metadata dictionary per PMID in ``pmids``.

    """
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


def text_or_none(node: ET.Element | None) -> str | None:
    """Return stripped text of an XML node if present."""
    if node is not None and node.text is not None:
        return node.text.strip()
    return None


def combine(items: Iterable[str]) -> str:
    """Combine non-empty items into a pipe-separated string without spaces."""
    return "|".join(x for x in items if x)


def find_one(node: ET.Element | None, xpath: str) -> ET.Element | None:
    """Safe wrapper around Element.find."""
    return node.find(xpath) if node is not None else None


def find_all(node: ET.Element | None, xpath: str) -> list[ET.Element]:
    """Safe wrapper around Element.findall returning a list."""
    return node.findall(xpath) if node is not None else []


def parse_pubmed_article(art: ET.Element) -> dict[str, Any]:
    """Parse PubMedArticle into a dictionary of selected fields."""
    mc = find_one(art, "./MedlineCitation")
    article = find_one(mc, "./Article") if mc is not None else None
    journal = find_one(article, "./Journal") if article is not None else None
    journal_issue = find_one(journal, "./JournalIssue") if journal is not None else None
    pagination = find_one(article, "./Pagination") if article is not None else None

    pmid = text_or_none(find_one(mc, "./PMID")) if mc is not None else None

    # DOI: prefer ArticleIdList[@IdType='doi'], fallback to ELocationID[@EIdType='doi'], normalize.
    doi = None
    if article is not None:
        for a in find_all(article, "./ArticleIdList/ArticleId"):
            if a.get("IdType") == "doi" and a.text:
                doi = a.text.strip()
                break
        if not doi:
            for el in find_all(article, "./ELocationID"):
                if el.get("EIdType") == "doi" and el.text:
                    doi = el.text.strip()
                    break
        if doi:
            lower = doi.lower()
            if lower.startswith("doi:"):
                doi = doi[4:].strip()
            for pref in ("https://doi.org/", "http://doi.org/", "doi.org/"):
                if doi.lower().startswith(pref):
                    doi = doi[len(pref) :].strip()
                    break

    article_title = (
        text_or_none(find_one(article, "./ArticleTitle"))
        if article is not None
        else None
    )

    # Abstract: join all segments, preserve label when present.
    article_abstract = None
    if article is not None:
        segments = find_all(article, "./Abstract/AbstractText")
        if segments:
            parts: list[str] = []
            for seg in segments:
                seg_text = text_or_none(seg)
                if seg_text:
                    label = seg.get("Label")
                    parts.append(f"{label}: {seg_text}" if label else seg_text)
            article_abstract = " ".join(parts) if parts else None
        if article_abstract is None:
            article_abstract = text_or_none(
                find_one(article, "./Abstract/AbstractText")
            )

    journal_title = text_or_none(find_one(journal, "./Title"))
    issn = text_or_none(find_one(journal, "./ISSN"))
    volume = text_or_none(find_one(journal_issue, "./Volume"))
    issue = text_or_none(find_one(journal_issue, "./Issue"))

    start_page = text_or_none(find_one(pagination, "./StartPage"))
    end_page = text_or_none(find_one(pagination, "./EndPage"))

    pubtypes = [
        p
        for p in (
            text_or_none(x)
            for x in find_all(article, "./PublicationTypeList/PublicationType")
        )
        if p
    ]

    mh_list = find_one(mc, "./MeshHeadingList")
    mesh_descriptors: list[str] = []
    mesh_qualifiers: list[str] = []
    if mh_list is not None:
        for mh in mh_list.findall("./MeshHeading"):
            d = text_or_none(find_one(mh, "./DescriptorName"))
            if d:
                mesh_descriptors.append(d)
            for q in mh.findall("./QualifierName"):
                qt = text_or_none(q)
                if qt:
                    mesh_qualifiers.append(qt)

    chemical_list: list[str] = []
    chem_list_node = find_one(mc, "./ChemicalList")
    if chem_list_node is not None:
        for chem in chem_list_node.findall("./Chemical"):
            nos = text_or_none(find_one(chem, "./NameOfSubstance"))
            if nos:
                chemical_list.append(nos)

    dr = find_one(mc, "./DateRevised")
    year_revised = text_or_none(find_one(dr, "./Year")) if dr is not None else ""
    month_revised = text_or_none(find_one(dr, "./Month")) if dr is not None else ""
    day_revised = text_or_none(find_one(dr, "./Day")) if dr is not None else ""

    dc = find_one(mc, "./DateCompleted")
    year_completed = text_or_none(find_one(dc, "./Year")) if dc is not None else ""
    month_completed = text_or_none(find_one(dc, "./Month")) if dc is not None else ""
    day_completed = text_or_none(find_one(dc, "./Day")) if dc is not None else ""

    return {
        "PubMed.PMID": pmid or "",
        "PubMed.DOI": doi or "",
        "PubMed.ArticleTitle": article_title or "-",
        "PubMed.Abstract": article_abstract or "-",
        "PubMed.JournalTitle": journal_title or "-",
        "PubMed.ISSN": issn or "",
        "PubMed.Volume": volume or "0",
        "PubMed.Issue": issue or "0",
        "PubMed.StartPage": start_page or "0",
        "PubMed.EndPage": end_page or "0",
        "PubMed.PublicationType": combine(pubtypes) or "unknown",
        "PubMed.MeSH_Descriptors": combine(mesh_descriptors) or "unknown",
        "PubMed.MeSH_Qualifiers": combine(mesh_qualifiers) or "unknown",
        "PubMed.ChemicalList": combine(chemical_list) or "unknown",
        "PubMed.DayRevised": day_revised or "0",
        "PubMed.MonthRevised": month_revised or "0",
        "PubMed.YearRevised": year_revised or "0",
        "PubMed.YearCompleted": year_completed or "0",
        "PubMed.MonthCompleted": month_completed or "0",
        "PubMed.DayCompleted": day_completed or "0",
    }


EMPTY_PUBMED: dict[str, str] = {
    "PubMed.PMID": "",
    "PubMed.DOI": "",
    "PubMed.ArticleTitle": "",
    "PubMed.Abstract": "",
    "PubMed.JournalTitle": "",
    "PubMed.ISSN": "",
    "PubMed.Volume": "",
    "PubMed.Issue": "",
    "PubMed.StartPage": "",
    "PubMed.EndPage": "",
    "PubMed.PublicationType": "",
    "PubMed.MeSH_Descriptors": "",
    "PubMed.MeSH_Qualifiers": "",
    "PubMed.ChemicalList": "",
    "PubMed.DayRevised": "",
    "PubMed.MonthRevised": "",
    "PubMed.YearRevised": "",
    "PubMed.YearCompleted": "",
    "PubMed.MonthCompleted": "",
    "PubMed.DayCompleted": "",
    "PubMed.Error": "",
}


def fetch_pubmed(
    session: requests.Session, pmid: str, sleep: float, cfg: PubMedCfg | None = None
) -> dict[str, str]:
    """Fetch metadata for a PMID from the PubMed API.

    Parameters
    ----------
    session:
        Active :class:`requests.Session`.
    pmid:
        PubMed identifier to query.
    sleep:
        Seconds to pause before making the request.
    cfg:
        Optional :class:`PubMedCfg` with API settings.
    """
    cfg = cfg or PubMedCfg()
    base = cfg.base.rstrip("/")
    url = f"{base}/efetch.fcgi?db=pubmed&id={pmid}&retmode=xml"
    timeout = (cfg.timeout_connect, cfg.timeout_read)
    text, error = _do_request(
        session, url, sleep, expect_json=False, retries=cfg.retries, timeout=timeout
    )
    result = EMPTY_PUBMED.copy()
    if error:
        result["PubMed.Error"] = error
        return result
    try:
        root = ET.fromstring(text)  # type: ignore[arg-type]
    except ET.ParseError as exc:
        result["PubMed.Error"] = str(exc)
        return result
    article = root.find(".//PubmedArticle")
    if article is None:
        result["PubMed.Error"] = "No PubmedArticle"
        return result
    result.update(parse_pubmed_article(article))
    return result


def fetch_semantic_scholar(
    session: requests.Session,
    pmid: str,
    sleep: float,
    cfg: SemanticScholarCfg | None = None,
) -> dict[str, str]:
    """Retrieve Semantic Scholar metadata for a single PMID.

    Parameters
    ----------
    session:
        Active :class:`requests.Session`.
    pmid:
        PubMed identifier to query.
    sleep:
        Seconds to pause before making the request.
    cfg:
        Optional :class:`SemanticScholarCfg` with API settings.

    Returns
    -------
    dict
        Mapping of Semantic Scholar fields and any error encountered.

    """
    fields = "publicationTypes,externalIds,paperId,venue"
    headers = {"Accept": "application/json"}
    cfg = cfg or SemanticScholarCfg()
    base = cfg.base.rstrip("/")
    url = f"{base}/paper/PMID:{pmid}"
    timeout = (cfg.timeout_connect, cfg.timeout_read)
    data, error = _do_request(
        session,
        url,
        sleep * 5,
        headers=headers,
        params={"fields": fields},
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
    pubtypes = data.get("publicationTypes") or []
    external_ids = data.get("externalIds") or {}
    doi = external_ids.get("DOI") or ""
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
    """Fetch metadata for multiple PMIDs using Semantic Scholar's batch API.

    Parameters
    ----------
    session:
        Active :class:`requests.Session`.
    pmids:
        List of PubMed identifiers.
    sleep:
        Seconds to pause before the request.
    cfg:
        Optional :class:`SemanticScholarCfg` with API settings.
    """
    if not pmids:
        return []

    cfg = cfg or SemanticScholarCfg()
    fields = "publicationTypes,externalIds,paperId,venue"
    headers = {"Accept": "application/json"}
    base = cfg.base.rstrip("/")
    url = f"{base}/paper/batch"

    prefixed_pmids = [f"PMID:{pmid}" for pmid in pmids]
    timeout = (cfg.timeout_connect, cfg.timeout_read)

    data, error = _do_request(
        session,
        url,
        sleep * 5,
        headers=headers,
        params={"fields": fields},
        json={"ids": prefixed_pmids},
        method="POST",
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

    results: list[dict[str, str]] = []
    if isinstance(data, list) and len(data) == len(pmids):
        for pmid, item in zip(pmids, data, strict=False):
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
    else:
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

    return results


def fetch_openalex(
    session: requests.Session,
    pmid: str,
    *,
    cfg: OpenAlexCfg,
    limiter: RateLimiter | None = None,
) -> dict[str, str]:
    """Retrieve OpenAlex metadata for ``pmid``.

    Parameters
    ----------
    session:
        Active :class:`requests.Session`.
    pmid:
        PubMed identifier to query.
    cfg:
        OpenAlex configuration providing base URL, timeouts and rate limits.
    limiter:
        Shared :class:`RateLimiter` instance. If ``None``, a limiter is
        retrieved via :func:`get_limiter` using the configuration values.


    Returns
    -------
    dict
        Mapping of OpenAlex fields and any error encountered.

    """

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
    """Retrieve Crossref metadata for a given DOI.

    Parameters
    ----------
    session:
        Active :class:`requests.Session`.
    doi:
        Digital Object Identifier to query.
    cfg:
        CrossRef configuration providing base URL, timeouts and rate limits.
    limiter:
        Shared :class:`RateLimiter` instance. If ``None``, a limiter is
        retrieved via :func:`get_limiter` using the configuration values.


    Returns
    -------
    dict
        Mapping of Crossref fields and any error encountered.

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


def print_results(records: list[dict[str, str]]) -> None:
    """Log result records instead of printing to ``stdout``.

    Parameters
    ----------
    records:
        Sequence of result dictionaries produced by ``main``.
    """
    log = logging.getLogger(__name__)
    try:
        from tabulate import tabulate

        use_tabulate = True
    except Exception:
        use_tabulate = False

    display_records = []
    for rec in records:
        d = rec.copy()
        title = (
            d.get("PubMed.ArticleTitle")
            or d.get("crossref.Title")
            or d.get("Title")
            or ""
        )
        d["Title"] = title[:77] + "..." if len(title) > 80 else title
        display_records.append(d)

    if use_tabulate:
        output = tabulate(display_records, headers="keys")
    else:
        output = json.dumps(display_records, ensure_ascii=False, indent=2)

    try:
        log.info(output)
    except UnicodeEncodeError:
        encoded = output.encode(sys.stdout.encoding or "utf-8", errors="replace")
        sys.stdout.buffer.write(encoded + b"\n")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Fetch publication metadata by PMID")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    parser.add_argument(
        "--input",
        dest="input_csv",
        default="input.csv",
        help="Input CSV path with PMID column",
    )
    parser.add_argument(
        "--output",
        dest="output_csv",
        default=None,
        help="Output CSV path (default: auto-generated)",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line interface entry point.

    Parameters
    ----------
    argv:
        Optional sequence of command-line arguments.

    Returns
    -------
    int
        Zero on success.
    """
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))
    log = logging.getLogger(__name__)

    cfg = Config()
    pubmed_cfg = cfg.pubmed
    semsch_cfg = cfg.semantic_scholar
    limiter = get_limiter("global", cfg.rate.global_rps, cfg.rate.global_burst)
    delay = 1.0 / cfg.rate.global_rps if cfg.rate.global_rps > 0 else 0.0

    cfg = Config()
    pubmed_cfg = cfg.pubmed
    semsch_cfg = cfg.semantic_scholar
    pmids = read_pmids(args.input, cfg=pubmed_cfg)
    openalex_limiter = get_limiter("openalex", cfg.openalex.rps, cfg.openalex.burst)
    crossref_limiter = get_limiter("crossref", cfg.crossref.rps, cfg.crossref.burst)

    records: list[dict[str, str]] = []
    batch_size = 100
    with session_with_retry(cfg.api, cfg.retry) as session:
        for i in range(0, len(pmids), batch_size):
            batch_pmids = pmids[i : i + batch_size]
            limiter.acquire()
            pubmed_list = fetch_pubmed_batch(
                session, batch_pmids, delay, cfg=pubmed_cfg
            )
            limiter.acquire()
            semsch_list = fetch_semantic_scholar_batch(
                session, batch_pmids, delay, cfg=semsch_cfg
            )
            semsch_map = {s.get("scholar.PMID"): s for s in semsch_list}
            for pubmed in pubmed_list:
                pmid = pubmed.get("PubMed.PMID", "")
                semsch = semsch_map.get(pmid, {})

                # Still fetching these individually

                openalex = fetch_openalex(
                    session, pmid, cfg=cfg.openalex, limiter=openalex_limiter
                )
                doi = pubmed.get("PubMed.DOI") or semsch.get("scholar.DOI") or ""
                crossref = fetch_crossref(
                    session, doi, cfg=cfg.crossref, limiter=crossref_limiter
                )

                combined: dict[str, str] = {}
                combined.update(pubmed)
                combined.update(semsch)
                combined.update(openalex)
                combined.update(crossref)
                print_results([combined])
                records.append(combined)

    df = pd.DataFrame.from_records(records)
    output_path = (
        Path(args.output_csv)
        if args.output_csv
        else Path(f"output_{Path(args.input_csv).stem}_{date.today():%Y%m%d}.csv")
    )
    write_csv_deterministic(df, output_path)
    log.info("written %s", output_path)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
