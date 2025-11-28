"""Crossref API client."""

from __future__ import annotations

from typing import Any
from urllib.parse import quote

import requests

from ..common.rate_limiter import RateLimiter
from ..config.models import CrossRefCfg, RetryCfg
from .base import RateLimitedRequestMixin

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

    base = cfg.base.rstrip("/")
    url = f"{base}/works/{quote(doi, safe='')}?mailto={quote(cfg.mailto)}"
    return RateLimitedRequestMixin.perform_request(
        session,
        url,
        cfg,
        limiter=limiter,
        retry_cfg=retry_cfg,
    )
