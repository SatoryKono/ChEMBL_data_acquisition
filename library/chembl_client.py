"""Shared HTTP utilities for ChEMBL API access."""

from __future__ import annotations

import random
import threading
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from typing import Any, cast

import requests
from cachetools import TTLCache  # type: ignore[import-untyped]
from requests import Session

from .config import ApiCfg, ChemblCfg, RetryCfg, session_with_retry
from .log import logger
from .rate_limiter import get_limiter, sleep


@dataclass
class ChemblClient:
    """HTTP client for the ChEMBL API with a TTL cache.

    Parameters
    ----------
    api:
        Global API settings providing the ``User-Agent`` header.
    retry:
        Retry configuration applied to all requests.
    chembl:
        Optional ChEMBL-specific configuration controlling cache TTL.
    session:
        Optional pre-configured :class:`requests.Session` instance; primarily
        intended for tests.
    """

    session: Session = field(init=False)
    cache: TTLCache[str, dict[str, Any]] = field(init=False)
    _cache_lock: threading.Lock = field(default_factory=threading.Lock, init=False)

    def __init__(
        self,
        api: ApiCfg | None = None,
        retry: RetryCfg | None = None,
        chembl: ChemblCfg | None = None,
        *,
        session: Session | None = None,
    ) -> None:
        api = api or ApiCfg()
        retry = retry or RetryCfg()
        self.session = session or session_with_retry(api, retry)
        ttl = chembl.cache_ttl if chembl is not None else ChemblCfg().cache_ttl
        self.cache = TTLCache(maxsize=1024, ttl=ttl)
        self._cache_lock = threading.Lock()

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
            If the HTTP request fails.
        ValueError
            If the response body is not valid JSON.
        """

        limiter = get_limiter("chembl", cfg.rps, cfg.burst)
        read_timeout = timeout if timeout is not None else cfg.timeout_read
        cache_key = url
        with self._cache_lock:
            cached = self.cache.get(cache_key)
            if cached is not None:
                logger.info("cache_hit", extra={"stage": "cache_hit", "url": url})
                return cast(dict[str, Any], cached)
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
                            return cast(dict[str, Any], cached)
                        self.cache[cache_key] = data
                        logger.info(
                            "cache_set", extra={"stage": "cache_set", "url": url}
                        )
                        return data
            except (requests.RequestException, ValueError) as exc:
                last_exc = exc
                if attempt >= cfg.retries:
                    logger.exception(
                        "request_fail", extra={"stage": "request_fail", "url": url}
                    )
                    break
                delay = cfg.backoff_factor * (2 ** (attempt - 1))
                delay += random.uniform(0, cfg.backoff_factor)
                sleep(delay)

        assert last_exc is not None
        raise last_exc

    def clear_cache(self) -> None:
        """Remove all entries from the in-memory cache."""

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
