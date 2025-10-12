"""Crossref API client."""

from __future__ import annotations

from typing import Any
from urllib.parse import quote

import requests

from ..common.rate_limiter import RateLimiter, get_limiter
from ..config.models import CrossRefCfg, RetryCfg
from .pubmed import _do_request

__all__ = ["fetch_crossref"]


def fetch_crossref(
    session: requests.Session,
    doi: str,
    *,
    cfg: CrossRefCfg,
    limiter: RateLimiter | None = None,
    retry_cfg: RetryCfg | None = None,
) -> tuple[dict[str, Any] | str | None, str]:
    """Request Crossref metadata for ``doi``."""

    if limiter is None:
        limiter = get_limiter("crossref", cfg.rps, cfg.burst)
    limiter.acquire()

    delay = 1 / cfg.rps if cfg.rps else 0
    base = cfg.base.rstrip("/")
    url = f"{base}/works/{quote(doi, safe='')}?mailto={quote(cfg.mailto)}"
    timeout = (cfg.timeout_connect, cfg.timeout_read)

    data, error = _do_request(
        session,
        url,
        delay,
        retries=cfg.retries,
        timeout=timeout,
        retry_cfg=retry_cfg,
    )
    return data, error
