"""Shared client helpers and mixins."""
from __future__ import annotations

from typing import Any

import requests

from ..common.log import logger
from ..common.rate_limiter import RateLimiter, get_limiter
from ..config.models import RetryCfg
from .pubmed import _do_request

__all__ = ["RateLimitedRequestMixin"]


class RateLimitedRequestMixin:
    """Utility methods for clients that need throttled HTTP requests."""

    @staticmethod
    def _resolve_delay(cfg: Any) -> float:
        rps = getattr(cfg, "rps", None)
        if rps:
            return 1.0 / rps
        delay = getattr(cfg, "delay", 0) or 0
        return float(delay)

    @staticmethod
    def _resolve_limiter(cfg: Any, limiter: RateLimiter | None) -> RateLimiter | None:
        if limiter is not None:
            return limiter
        rps = getattr(cfg, "rps", None)
        if rps is None or rps <= 0:
            return None
        burst = getattr(cfg, "burst", None)
        name = getattr(cfg, "limiter_name", None) or getattr(cfg, "base", None)
        if name is None:
            name = cfg.__class__.__name__.lower()
        return get_limiter(str(name), rps, burst)

    @staticmethod
    def perform_request(
        session: requests.Session,
        url: str,
        cfg: Any,
        *,
        limiter: RateLimiter | None = None,
        retry_cfg: RetryCfg | None = None,
        headers: dict[str, str] | None = None,
        params: dict[str, Any] | None = None,
        json: Any = None,
        method: str = "GET",
    ) -> tuple[dict[str, Any] | str | None, str]:
        effective_limiter = RateLimitedRequestMixin._resolve_limiter(cfg, limiter)
        if effective_limiter is not None:
            effective_limiter.acquire()

        delay = RateLimitedRequestMixin._resolve_delay(cfg)
        timeout = (cfg.timeout_connect, cfg.timeout_read)

        data, error = _do_request(
            session,
            url,
            delay,
            headers=headers,
            params=params,
            json=json,
            method=method,
            retries=cfg.retries,
            timeout=timeout,
            retry_cfg=retry_cfg,
        )
        if error:
            logger.debug("request_fail", extra={"stage": "request_fail", "url": url})
        else:
            logger.debug("request_ok", extra={"stage": "request_ok", "url": url})
        return data, error
