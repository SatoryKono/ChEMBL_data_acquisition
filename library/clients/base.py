"""Reusable HTTP client helpers for external data providers."""

from __future__ import annotations

import random
import threading
from dataclasses import dataclass, field
from types import TracebackType
from typing import Any

import requests
from cachetools import TTLCache
from requests import Session

from ..config import ApiCfg, ChemblCacheCfg, RetryCfg, session_with_retry
from ..rate_limiter import get_limiter, sleep
from ..utils.logging import logger


@dataclass
class ChemblClient:
  """HTTP client for the ChEMBL API with a TTL cache."""

  session: Session = field(init=False)
  cache: TTLCache[str, dict[str, Any]] = field(init=False)
  _cache_lock: threading.Lock = field(default_factory=threading.Lock, init=False)

  def __init__(
    self,
    api: ApiCfg | None = None,
    retry: RetryCfg | None = None,
    chembl: ChemblCacheCfg | None = None,
    *,
    session: Session | None = None,
  ) -> None:
    api = api or ApiCfg(user_agent="chembl-da/0.1 (mailto:contact@example.org)")
    retry = retry or RetryCfg()
    self.session = session or session_with_retry(api, retry)
    ttl = chembl.cache_ttl if chembl is not None else ChemblCacheCfg().cache_ttl
    maxsize = (
      chembl.cache_maxsize
      if chembl is not None
      else ChemblCacheCfg().cache_maxsize
    )
    self.cache = TTLCache(maxsize=maxsize, ttl=ttl)
    self._cache_lock = threading.Lock()

  def close(self) -> None:
    """Close the underlying HTTP session."""

    self.session.close()

  def __enter__(self) -> "ChemblClient":
    return self

  def __exit__(
    self,
    exc_type: type[BaseException] | None,
    exc: BaseException | None,
    tb: TracebackType | None,
  ) -> None:
    self.close()

  def request_json(
    self,
    url: str,
    *,
    cfg: ApiCfg,
    timeout: float | None = None,
  ) -> dict[str, Any]:
    """Return JSON content from ``url`` using retry and rate limiting."""

    limiter = get_limiter("chembl", cfg.rps, cfg.burst)
    read_timeout = timeout if timeout is not None else cfg.timeout_read
    cache_key = url
    with self._cache_lock:
      cached = self.cache.get(cache_key)
      if cached is not None:
        logger.info(
          "cache_hit", extra={"url": url, "rps": cfg.rps, "status": "hit"}
        )
        return cached
      logger.info(
        "cache_miss", extra={"url": url, "rps": cfg.rps, "status": "miss"}
      )

    last_exc: requests.RequestException | ValueError | None = None
    total_attempts = cfg.retries + 1

    for attempt in range(1, total_attempts + 1):
      limiter.acquire()
      event = "request_start" if attempt == 1 else "request_retry"
      logger.info(event, extra={"url": url, "attempt": attempt, "rps": cfg.rps})
      try:
        with self.session.get(
          url, timeout=(cfg.timeout_connect, read_timeout)
        ) as response:
          response.raise_for_status()
          try:
            data = response.json()
          except ValueError as exc:
            logger.exception("json_error", extra={"url": url})
            raise ValueError(f"invalid JSON in response from {url}") from exc
          logger.info(
            "request_ok",
            extra={
              "url": url,
              "status": getattr(response, "status_code", None),
              "rps": cfg.rps,
            },
          )
          with self._cache_lock:
            cached = self.cache.get(cache_key)
            if cached is not None:
              return cached
            self.cache[cache_key] = data
            logger.info("cache_set", extra={"url": url, "rps": cfg.rps})
            return data
      except (requests.RequestException, ValueError) as exc:
        last_exc = exc
        if attempt >= total_attempts:
          logger.exception(
            "request_fail",
            extra={"url": url, "status": None, "rps": cfg.rps},
          )
          break
        delay = cfg.backoff_factor * (2 ** (attempt - 1))
        delay += random.uniform(0, cfg.backoff_factor)
        sleep(delay)

    if last_exc is not None:
      raise last_exc
    raise RuntimeError(f"Request loop exited unexpectedly for {url}")

  def clear_cache(self) -> None:
    """Remove all entries from the in-memory cache."""

    with self._cache_lock:
      self.cache.clear()


__all__ = ["ChemblClient"]
