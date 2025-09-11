"""Shared HTTP utilities for ChEMBL API access."""

from __future__ import annotations

from typing import Any, Dict, Iterable, Iterator

import random
import requests
from requests import Session

from .config import ApiCfg, RetryCfg, session_with_retry
from .rate_limiter import sleep
from .log import logger

_CACHE: Dict[str, dict[str, Any]] = {}

_session: Session = session_with_retry(ApiCfg(), RetryCfg())


def init_session(api: ApiCfg, retry: RetryCfg) -> None:
    """Initialise the shared HTTP session.

    Parameters
    ----------
    api:
        Global API settings providing the ``User-Agent`` header.
    retry:
        Retry configuration applied to all requests.
    """

    global _session
    _session = session_with_retry(api, retry)


def request_json(
    url: str, *, cfg: ApiCfg, timeout: float | None = None
) -> dict[str, Any]:
    """Return JSON content from *url*.

    Parameters
    ----------
    url:
        API endpoint to query.
    cfg:
        Configuration providing timeout and retry settings.
    timeout:
        Optional override for the read timeout in seconds.

    Returns
    -------
    dict[str, Any]
        Parsed JSON document.

    Raises
    ------
    requests.RequestException
        If the HTTP request fails.
    ValueError
        If the response body is not valid JSON.

    """
    read_timeout = timeout if timeout is not None else cfg.timeout_read
    cache_key = url
    if cache_key in _CACHE:
        logger.info("cache_hit", extra={"stage": "cache_hit", "url": url})
        return _CACHE[cache_key]
    logger.info("cache_miss", extra={"stage": "cache_miss", "url": url})

    for attempt in range(1, cfg.retries + 1):
        event = "request_start" if attempt == 1 else "request_retry"
        logger.info(event, extra={"stage": event, "url": url, "attempt": attempt})
        try:
            with _session.get(
                url, timeout=(cfg.timeout_connect, read_timeout)
            ) as response:
                response.raise_for_status()
                data = response.json()
                logger.info(
                    "request_ok",
                    extra={
                        "stage": "request_ok",
                        "url": url,
                        "status": getattr(response, "status_code", None),
                    },
                )
                _CACHE[cache_key] = data
                logger.info("cache_set", extra={"stage": "cache_set", "url": url})
                return data
        except (requests.RequestException, ValueError):
            if attempt >= cfg.retries:
                logger.exception(
                    "request_fail", extra={"stage": "request_fail", "url": url}
                )
                raise
            # Exponential backoff with jitter to avoid thundering herd problems
            delay = cfg.backoff_factor * (2 ** (attempt - 1))
            delay += random.uniform(0, cfg.backoff_factor)
            sleep(delay)
    raise requests.RequestException(f"request_json failed for url: {url}")


def _chunked(items: Iterable[str], size: int) -> Iterator[list[str]]:
    """Yield ``size``-sized lists from *items*.

    Parameters
    ----------
    items:
        Iterable of identifiers to split.
    size:
        Desired chunk size; must be positive.

    Yields
    ------
    list[str]
        Subsequences of ``items`` with at most ``size`` elements.

    Raises
    ------
    ValueError
        If ``size`` is not a positive integer.

    """
    if size <= 0:
        raise ValueError("size must be a positive integer")

    chunk: list[str] = []
    for item in items:
        chunk.append(item)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk
