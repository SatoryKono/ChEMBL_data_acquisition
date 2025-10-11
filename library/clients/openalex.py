"""OpenAlex API client."""

from __future__ import annotations

from typing import Any
from urllib.parse import quote

import requests

from ..common.log import logger
from ..common.rate_limiter import RateLimiter, get_limiter
from ..config.models import OpenAlexCfg, RetryCfg
from .pubmed import _do_request

__all__ = ["fetch_openalex"]


def fetch_openalex(
    session: requests.Session,
    pmid: str,
    *,
    cfg: OpenAlexCfg,
    limiter: RateLimiter | None = None,
    retry_cfg: RetryCfg | None = None,
) -> tuple[dict[str, Any] | str | None, str]:
    """Request OpenAlex metadata for ``pmid``."""

    if limiter is None:
        limiter = get_limiter("openalex", cfg.rps, cfg.burst)
    limiter.acquire()

    delay = 1 / cfg.rps if cfg.rps else 0
    base = cfg.base.rstrip("/")
    url = f"{base}/works/pmid:{pmid}?mailto={quote(cfg.mailto)}"
    timeout = (cfg.timeout_connect, cfg.timeout_read)

    data, error = _do_request(
        session,
        url,
        delay,
        retries=cfg.retries,
        timeout=timeout,
        retry_cfg=retry_cfg,
    )
    if error:
        logger.debug("request_fail", extra={"stage": "request_fail", "url": url})
    else:
        logger.debug("request_ok", extra={"stage": "request_ok", "url": url})
    return data, error
