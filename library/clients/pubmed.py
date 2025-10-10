"""Low-level HTTP client for PubMed E-Utilities."""

from __future__ import annotations


import random

from collections.abc import Callable, Mapping
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime

from typing import Any

import requests

from ..common.fetch_retry import compute_backoff_delay
from ..config.models import PubMedCfg, RetryCfg
from ..common.log import logger
from ..common.rate_limiter import sleep

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
) -> tuple[int, str, dict[str, Any] | str | None, str, Mapping[str, str]]:
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
                return status_code, text, None, str(exc), resp.headers
            return status_code, text, content, "", resp.headers
        return status_code, text, text or "", "", resp.headers


def _retry_after_delay(headers: Mapping[str, str]) -> float | None:
    """Parse the ``Retry-After`` header into a delay in seconds."""

    header = headers.get("Retry-After")
    if header is None:
        return None
    value = header.strip()
    if not value:
        return None
    try:
        seconds = float(value)
    except ValueError:
        try:
            retry_at = parsedate_to_datetime(value)
        except (TypeError, ValueError, OverflowError):
            return None
        if retry_at.tzinfo is None:
            retry_at = retry_at.replace(tzinfo=timezone.utc)
        now = datetime.now(timezone.utc)
        return max((retry_at - now).total_seconds(), 0.0)
    else:
        return max(seconds, 0.0)


def _handle_response(
    url: str,
    status_code: int,
    text: str,
    content: dict[str, Any] | str | None,
    parse_error: str,
    expect_json: bool,
    attempt: int,
    retries: int,
    headers: Mapping[str, str],
) -> tuple[dict[str, Any] | str | None, str, bool, float | None]:
    """Process an HTTP response and decide whether to retry."""

    if status_code in (429, 500, 502, 503, 504):
        retry_after = _retry_after_delay(headers)
        if attempt >= retries:
            logger.info(
                "request_fail",
                extra={"stage": "request_fail", "url": url, "status": status_code},
            )
            return None, f"HTTP {status_code}: {text[:100]}", False, retry_after
        return None, "", True, retry_after
    if status_code == 404:
        logger.info(
            "request_fail",
            extra={"stage": "request_fail", "url": url, "status": status_code},
        )
        return None, "PMID not found", False, None
    if status_code == 400:
        logger.info(
            "request_fail",
            extra={"stage": "request_fail", "url": url, "status": status_code},
        )
        return None, f"Bad request: {text[:100]}", False, None
    if status_code != 200:
        logger.info(
            "request_fail",
            extra={"stage": "request_fail", "url": url, "status": status_code},
        )
        return None, f"HTTP {status_code}: {text[:100]}", False, None
    if expect_json:
        if parse_error:
            logger.info("request_fail", extra={"stage": "request_fail", "url": url})
            return None, f"Invalid JSON: {parse_error}", False, None
        logger.info(
            "request_ok",
            extra={"stage": "request_ok", "url": url, "status": status_code},
        )
        return content, "", False, None
    logger.info(
        "request_ok",
        extra={"stage": "request_ok", "url": url, "status": status_code},
    )
    return content or "", "", False, None


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
    *,
    jitter: Callable[[float], float] | None = None,
) -> float:
    """Calculate the delay before the next retry attempt."""

    if attempt <= 0:
        return 0.0

    min_delay = max(base_delay, 0.0)
    delay = min_delay

    if retry_cfg is not None and retry_cfg.backoff_factor > 0:
        backoff = compute_backoff_delay(attempt, retry_cfg, jitter=jitter)
        candidate = backoff
        if jitter is None:
            jitter_value = random.uniform(0.0, retry_cfg.backoff_factor)
            candidate = backoff + jitter_value
        if retry_cfg.backoff_cap is not None:
            candidate = min(candidate, retry_cfg.backoff_cap)
        delay = max(delay, candidate)

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
    jitter: Callable[[float], float] | None = None,
    **kwargs: Any,
) -> tuple[dict[str, Any] | str | None, str]:
    """Perform an HTTP request with retry logic."""

    active_retry_cfg = retry_cfg or RetryCfg()
    effective_jitter = jitter
    if effective_jitter is None:
        effective_jitter = active_retry_cfg.build_jitter()
    retry_after_delay: float | None = None
    for attempt in range(retries + 1):
        event = "request_start" if attempt == 0 else "request_retry"
        extra = {"stage": event, "url": url, "attempt": attempt + 1}

        if attempt:
            header_delay = retry_after_delay
            retry_after_delay = None
            retry_delay = (
                header_delay
                if header_delay is not None
                else _retry_delay(
                    attempt,
                    delay,
                    active_retry_cfg,
                    timeout,
                    jitter=effective_jitter,
                )
            )
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
            status_code, text, content, parse_error, headers = _make_request(
                session, url, expect_json, method, timeout, **kwargs
            )
        except requests.RequestException as exc:
            if attempt >= retries:  # pragma: no cover - network errors
                logger.warning(
                    "request_fail",
                    extra={
                        "stage": "request_fail",
                        "url": url,
                        "error": str(exc),
                    },
                )
                return None, str(exc)
            continue

        data, error, retry, retry_after_delay = _handle_response(
            url,
            status_code,
            text,
            content,
            parse_error,
            expect_json,
            attempt,
            retries,
            headers,
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
    jitter: Callable[[float], float] | None = None,
) -> tuple[str | None, str]:
    """Fetch raw PubMed XML for ``pmids`` using a single request."""

    cfg = cfg or PubMedCfg()
    ids = ",".join(pmids)
    base = cfg.base.rstrip("/")
    params: dict[str, str] = {
        "db": "pubmed",
        "id": ids,
        "retmode": "xml",
    }
    if cfg.tool:
        params["tool"] = cfg.tool
    if cfg.email:
        params["email"] = cfg.email
    url = f"{base}/efetch.fcgi"
    timeout = (cfg.timeout_connect, cfg.timeout_read)
    effective_jitter = jitter
    if effective_jitter is None and retry_cfg is not None:
        effective_jitter = retry_cfg.build_jitter()

    text, error = _do_request(
        session,
        url,
        sleep,
        expect_json=False,
        retries=cfg.retries,
        timeout=timeout,
        retry_cfg=retry_cfg,
        jitter=effective_jitter,
        params=params,
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
    jitter: Callable[[float], float] | None = None,
) -> tuple[str | None, str]:
    """Fetch raw PubMed XML for a single PMID."""

    text, error = fetch_pubmed_batch(
        session,
        [pmid],
        sleep,
        cfg=cfg,
        retry_cfg=retry_cfg,
        jitter=jitter,
    )
    if error:
        return None, error
    return text, ""


class PubMedClient:
    """Thin wrapper around the PubMed HTTP helpers."""

    def __init__(
        self,
        cfg: PubMedCfg | None = None,
        *,
        jitter: Callable[[float], float] | None = None,
    ) -> None:
        self.cfg = cfg or PubMedCfg()
        self._jitter = jitter

    def fetch_pubmed_batch(
        self,
        session: requests.Session,
        pmids: list[str],
        sleep: float,
        *,
        retry_cfg: RetryCfg | None = None,
    ) -> tuple[str | None, str]:
        """Retrieve raw XML for ``pmids`` using configured settings."""

        effective_jitter = self._jitter
        if effective_jitter is None and retry_cfg is not None:
            effective_jitter = retry_cfg.build_jitter()

        return fetch_pubmed_batch(
            session,
            pmids,
            sleep,
            cfg=self.cfg,
            retry_cfg=retry_cfg,
            jitter=effective_jitter,
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

        effective_jitter = self._jitter
        if effective_jitter is None and retry_cfg is not None:
            effective_jitter = retry_cfg.build_jitter()

        return fetch_pubmed(
            session,
            pmid,
            sleep,
            cfg=self.cfg,
            retry_cfg=retry_cfg,
            jitter=effective_jitter,
        )
