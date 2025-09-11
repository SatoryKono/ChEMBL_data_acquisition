"""Shared HTTP utilities for ChEMBL API access.

This module provides :class:`ChemblClient`, a lightweight wrapper around
:class:`requests.Session` equipped with an in-memory cache.  The previous
implementation relied on module level globals which caused state to leak
between tests and CLI invocations.  ``ChemblClient`` encapsulates the
necessary state so each consumer can manage its own session and cache.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Iterator, cast

import random
import threading
import requests
from requests import Session
from cachetools import TTLCache  # type: ignore[import-untyped]

from .config import ApiCfg, ChemblCfg, RetryCfg, session_with_retry
from .rate_limiter import get_limiter, sleep
from .log import logger


@dataclass
class ChemblClient:
    """HTTP client with caching for the ChEMBL API.

    Parameters
    ----------
    api:
        Global API settings providing headers and rate limits.  A default
        instance is created when omitted.
    retry:
        Retry configuration for the underlying HTTP session.  A default is
        used when omitted.
    chembl:
        Optional ChEMBL specific settings controlling cache TTL.
    session:
        Pre-configured :class:`requests.Session` instance.  Primarily useful
        for tests where a dummy session should be injected.
    cache:
        Optional external cache.  When ``None`` a :class:`cachetools.TTLCache`
        is created using ``chembl.cache_ttl``.
    """

    session: Session
    cache: TTLCache[str, dict[str, Any]]

    def __init__(
        self,
        api: ApiCfg | None = None,
        retry: RetryCfg | None = None,
        chembl: ChemblCfg | None = None,
        *,
        session: Session | None = None,
        cache: TTLCache[str, dict[str, Any]] | None = None,
    ) -> None:
        api_cfg = api or ApiCfg()
        retry_cfg = retry or RetryCfg()
        chembl_cfg = chembl or ChemblCfg()
        self.session = session or session_with_retry(api_cfg, retry_cfg)
        self.cache = cache or TTLCache(maxsize=1024, ttl=chembl_cfg.cache_ttl)
        # Locks guard access to ``cache`` as cachetools' TTLCache is not thread
        # safe by default.  Each instance manages its own lock to avoid shared
        # state between clients.
        self._cache_lock = threading.Lock()

    # ------------------------------------------------------------------
    # Request handling
    # ------------------------------------------------------------------
    def request_json(
        self, url: str, *, cfg: ApiCfg, timeout: float | None = None
    ) -> dict[str, Any]:
        """Return JSON content from ``url``.

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
            If the HTTP request fails after exhausting retries.
        ValueError
            If the response body cannot be decoded as JSON.
        """

        limiter = get_limiter("chembl", cfg.rps, cfg.burst)
        read_timeout = timeout if timeout is not None else cfg.timeout_read
        cache_key = url
        with self._cache_lock:
            cached = self.cache.get(cache_key)
            if cached is not None:
                logger.info("cache_hit", extra={"stage": "cache_hit", "url": url})
                return cached
            logger.info("cache_miss", extra={"stage": "cache_miss", "url": url})

        last_exc: requests.RequestException | ValueError | None = None
        for attempt in range(1, cfg.retries + 1):
            limiter.acquire()
            event = "request_start" if attempt == 1 else "request_retry"
            logger.info(event, extra={"stage": event, "url": url, "attempt": attempt})
            try:
                with self.session.get(
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
                    with self._cache_lock:
                        cached = self.cache.get(cache_key)
                        if cached is not None:
                            return cached
                        self.cache[cache_key] = data
                        logger.info(
                            "cache_set", extra={"stage": "cache_set", "url": url}
                        )
                    return data
            except (
                requests.RequestException,
                ValueError,
            ) as exc:  # pragma: no cover - network errors
                last_exc = exc
                if attempt >= cfg.retries:
                    logger.exception(
                        "request_fail", extra={"stage": "request_fail", "url": url}
                    )
                    break
                # Exponential backoff with jitter to avoid thundering herd issues
                delay = cfg.backoff_factor * (2 ** (attempt - 1))
                delay += random.uniform(0, cfg.backoff_factor)
                sleep(delay)

        assert last_exc is not None
        raise last_exc

    def clear_cache(self) -> None:
        """Remove all entries from the in-memory cache.

        This helper is primarily intended for tests to avoid interference
        from previously cached responses.
        """

        with self._cache_lock:
            self.cache.clear()


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


__all__ = ["ChemblClient", "_chunked"]
