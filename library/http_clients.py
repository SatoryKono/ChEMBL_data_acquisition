"""HTTP clients for external metadata services.

This module contains low-level helpers that query the OpenAlex and CrossRef
APIs.  Higher level modules such as :mod:`library.pubmed_library` and
:mod:`library.openalex_crossref_library` import these functions to avoid a
circular dependency between them.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any
from urllib.parse import quote

import requests

from .config import CrossRefCfg, OpenAlexCfg
from .log import logger
from .rate_limiter import RateLimiter, get_limiter, sleep


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
        Active :class:`requests.Session` used to issue the call.
    url:
        Endpoint to query.
    delay:
        Base number of seconds to wait between attempts.
    expect_json:
        Whether to parse the response as JSON.
    retries:
        Number of additional attempts after the initial one.
    method:
        HTTP method to use, either ``"GET"`` or ``"POST"``.
    timeout:
        Maximum seconds to wait for each HTTP request. May be a single float or
        a ``(connect, read)`` tuple.
    **kwargs:
        Additional parameters passed to ``session.get`` or ``session.post``.

    Returns
    -------
    tuple
        Pair of parsed data (``dict``/``str``/``None``) and an error message.
        The error string is empty on success.
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


def fetch_openalex(
    session: requests.Session,
    pmid: str,
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
    combine = "|".join
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
        Mapping of CrossRef fields and any error encountered.
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


__all__ = ["fetch_openalex", "fetch_crossref", "_do_request"]
