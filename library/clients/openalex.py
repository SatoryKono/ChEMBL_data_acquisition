"""OpenAlex API client."""

from __future__ import annotations

from typing import Any
from urllib.parse import quote

import requests

from ..common.rate_limiter import RateLimiter
from ..config.models import OpenAlexCfg, RetryCfg
from .base import RateLimitedRequestMixin

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

    base = cfg.base.rstrip("/")
    url = f"{base}/works/pmid:{pmid}?mailto={quote(cfg.mailto)}"
    return RateLimitedRequestMixin.perform_request(
        session,
        url,
        cfg,
        limiter=limiter,
        retry_cfg=retry_cfg,
    )
