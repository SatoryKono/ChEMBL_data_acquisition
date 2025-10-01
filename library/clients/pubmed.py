"""Low-level HTTP client for PubMed E-Utilities."""

from __future__ import annotations

import random
from collections.abc import Callable
from typing import Any

import requests

from ..config import PubMedCfg, RetryCfg
from ..log import logger
from ..rate_limiter import sleep

__all__ = [
    "PubMedClient",
    "_make_request",
    "_handle_response",
    "_do_request",
    "fetch_pubmed_batch",
    "fetch_pubmed",
]


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


def _max_timeout(timeout: float | tuple[float, float] | None) -> float | None:
    """Return the maximum timeout value from ``timeout``."""

    if timeout is None:
        return None
    if isinstance(timeout, tuple):
        values = [float(value) for value in timeout if value is not None]
        if not values:
            return None
        return max(values)
    return float(timeout)


def _retry_delay(
    attempt: int,
    base_delay: float,
    retry_cfg: RetryCfg | None,
    timeout: float | tuple[float, float] | None,
) -> float:
    """Calculate the delay before the next retry attempt."""

    if attempt <= 0:
        return 0.0

    min_delay = max(base_delay, 0.0)
    delay = min_delay

    if retry_cfg is not None and retry_cfg.backoff_factor > 0:
        backoff = retry_cfg.backoff_factor * (2 ** (attempt - 1))
        jitter = random.uniform(0.0, retry_cfg.backoff_factor)
        delay = max(delay, backoff + jitter)

    timeout_cap = _max_timeout(timeout)
    if timeout_cap is not None:
        delay = min(delay, timeout_cap)

    return delay


def _do_request(
    session: requests.Session,
    url: str,
    delay: float,
    expect_json: bool = True,
    retries: int = 2,
    method: str = "GET",
    timeout: float | tuple[float, float] = 10,
    retry_cfg: RetryCfg | None = None,
    **kwargs: Any,
) -> tuple[dict[str, Any] | str | None, str]:
    """Perform an HTTP request with retry logic."""

    for attempt in range(retries + 1):
        event = "request_start" if attempt == 0 else "request_retry"
        extra = {"stage": event, "url": url, "attempt": attempt + 1}

        if attempt:
            retry_delay = _retry_delay(attempt, delay, retry_cfg, timeout)
            extra["delay"] = retry_delay
            logger.info(event, extra=extra)
            if retry_delay > 0:
                logger.debug(
                    "retry_sleep",
                    extra={
                        "url": url,
                        "attempt": attempt + 1,
                        "delay": retry_delay,
                    },
                )
                sleep(retry_delay)
        else:
            logger.info(event, extra=extra)

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
    *,
    retry_cfg: RetryCfg | None = None,
) -> tuple[str | None, str]:
    """Fetch raw PubMed XML for ``pmids`` using a single request."""

    cfg = cfg or PubMedCfg()
    ids = ",".join(pmids)
    base = cfg.base.rstrip("/")
    url = f"{base}/efetch.fcgi?db=pubmed&id={ids}&retmode=xml"
    timeout = (cfg.timeout_connect, cfg.timeout_read)
    text, error = _do_request(
        session,
        url,
        sleep,
        expect_json=False,
        retries=cfg.retries,
        timeout=timeout,
        retry_cfg=retry_cfg,
    )
    if error:
        return None, error
    if not isinstance(text, str):
        return None, "Invalid PubMed response"
    return text, ""


def fetch_pubmed(
    session: requests.Session,
    pmid: str,
    sleep: float,
    cfg: PubMedCfg | None = None,
    *,
    retry_cfg: RetryCfg | None = None,
) -> tuple[str | None, str]:
    """Fetch raw PubMed XML for a single PMID."""

    text, error = fetch_pubmed_batch(
        session, [pmid], sleep, cfg=cfg, retry_cfg=retry_cfg
    )
    if error:
        return None, error
    return text, ""


class PubMedClient:
    """Thin wrapper around the PubMed HTTP helpers."""

    def __init__(self, cfg: PubMedCfg | None = None) -> None:
        self.cfg = cfg or PubMedCfg()

    def fetch_pubmed_batch(
        self,
        session: requests.Session,
        pmids: list[str],
        sleep: float,
        *,
        retry_cfg: RetryCfg | None = None,
    ) -> tuple[str | None, str]:
        """Retrieve raw XML for ``pmids`` using configured settings."""

        return fetch_pubmed_batch(
            session, pmids, sleep, cfg=self.cfg, retry_cfg=retry_cfg
        )

    def fetch_pubmed(
        self,
        session: requests.Session,
        pmid: str,
        sleep: float,
        *,
        retry_cfg: RetryCfg | None = None,
    ) -> tuple[str | None, str]:
        """Retrieve raw XML for a single PMID."""

        return fetch_pubmed(
            session, pmid, sleep, cfg=self.cfg, retry_cfg=retry_cfg
        )
