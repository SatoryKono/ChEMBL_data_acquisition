"""Shared HTTP utilities for ChEMBL API access."""

from __future__ import annotations


from typing import Any, Iterable, Iterator, cast


import random
import threading
import requests
from requests import Session

from cachetools import TTLCache  # type: ignore[import-untyped]

from .config import ApiCfg, RetryCfg, session_with_retry
from .rate_limiter import get_limiter, sleep
from .log import logger

# Cache entries expire after one hour to avoid serving stale data. The TTL can
# be adjusted in the future via configuration if required.
_CACHE: TTLCache[str, dict[str, Any]] = TTLCache(maxsize=1024, ttl=3600)

_session: Session | None = None
_session_lock = threading.Lock()


def init_session(api: ApiCfg, retry: RetryCfg) -> None:
    """Initialise the shared HTTP session.

    The provided ``api`` and ``retry`` configurations are forwarded to
    :func:`session_with_retry`, ensuring that subsequent requests use the
    correct ``User-Agent`` and retry policy.

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
    limiter = get_limiter("chembl", cfg.rps, cfg.burst)
    read_timeout = timeout if timeout is not None else cfg.timeout_read
    cache_key = url
    if cache_key in _CACHE:
        logger.info("cache_hit", extra={"stage": "cache_hit", "url": url})
        return _CACHE[cache_key]
    logger.info("cache_miss", extra={"stage": "cache_miss", "url": url})

    global _session
    if _session is None:
        with _session_lock:
            if _session is None:
                _session = session_with_retry(ApiCfg(), RetryCfg())
    session = _session
    assert session is not None  # noqa: S101 - ensure session exists

    last_exc: requests.RequestException | ValueError | None = None

    for attempt in range(1, cfg.retries + 1):
        limiter.acquire()
        event = "request_start" if attempt == 1 else "request_retry"
        logger.info(event, extra={"stage": event, "url": url, "attempt": attempt})
        try:
            with session.get(
                url, timeout=(cfg.timeout_connect, read_timeout)
            ) as response:
                response.raise_for_status()
                data: dict[str, Any] = cast(dict[str, Any], response.json())
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
        except (requests.RequestException, ValueError) as exc:
            last_exc = exc
            if attempt >= cfg.retries:
                logger.exception(
                    "request_fail", extra={"stage": "request_fail", "url": url}
                )
                break
            # Exponential backoff with jitter to avoid thundering herd problems
            delay = cfg.backoff_factor * (2 ** (attempt - 1))
            delay += random.uniform(0, cfg.backoff_factor)
            sleep(delay)

    assert last_exc is not None
    raise last_exc


def clear_cache() -> None:
    """Remove all entries from the in-memory cache.

    This helper is primarily intended for tests to avoid interference
    from previously cached responses.
    """

    _CACHE.clear()


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
