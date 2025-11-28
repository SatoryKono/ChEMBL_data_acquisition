from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from typing import Any

import requests

from ..common.log import logger
from ..common.rate_limiter import RateLimiter, get_limiter
from .pubmed import _do_request


@dataclass(slots=True)
class BasePaginatedClient(ABC):
    """Base class for clients supporting limit/offset pagination."""

    limit_param: str = "limit"
    offset_param: str = "offset"
    page_size: int = 100

    def __post_init__(self) -> None:
        if self.page_size <= 0:
            raise ValueError("page_size must be positive")

    def _build_pagination_params(
        self,
        *,
        params: Mapping[str, Any] | None,
        limit: int,
        offset: int,
    ) -> dict[str, Any]:
        if limit <= 0:
            raise ValueError("limit must be positive")
        if offset < 0:
            raise ValueError("offset must be non-negative")

        effective_params: dict[str, Any] = {
            self.limit_param: limit,
            self.offset_param: offset,
        }
        if params:
            effective_params.update(params)
        return effective_params

    def iter_pages(
        self,
        endpoint: str,
        *,
        params: Mapping[str, Any] | None = None,
        start_offset: int = 0,
    ) -> Iterator[Mapping[str, Any]]:
        """Yield subsequent pages until the API stops providing a next link."""

        offset = start_offset
        if offset < 0:
            raise ValueError("start_offset must be non-negative")

        while True:
            page_params = self._build_pagination_params(
                params=params, limit=self.page_size, offset=offset
            )
            payload = self.request_json(endpoint, params=page_params)
            yield payload

            page_meta = payload.get("page_meta") if isinstance(payload, Mapping) else None
            next_url = page_meta.get("next") if isinstance(page_meta, Mapping) else None
            if not next_url:
                break
            offset += self.page_size

    @abstractmethod
    def request_json(
        self,
        endpoint: str,
        *,
        params: Mapping[str, Any] | None = None,
        **kwargs: Any,
    ) -> Mapping[str, Any]:
        """Perform a JSON request and return a response payload."""


class RateLimitedRequestMixin:
    """Provide a reusable rate-limited HTTP request helper."""

    limiter_name: str | None = None

    def __init__(self, limiter_name: str | None = None) -> None:
        if limiter_name is not None:
            self.limiter_name = limiter_name

    def _resolve_limiter_name(self, cfg: Any) -> str:
        if self.limiter_name:
            return self.limiter_name
        return cfg.__class__.__name__.lower()

    def _resolve_delay(self, cfg: Any) -> float:
        rps = getattr(cfg, "rps", None)
        if rps:
            return 1.0 / rps
        return float(getattr(cfg, "delay", 0.0) or 0.0)

    def perform_request(
        self,
        session: requests.Session,
        url: str,
        cfg: Any,
        *,
        limiter: RateLimiter | None = None,
        retry_cfg: Any | None = None,
        headers: Mapping[str, str] | None = None,
        params: Mapping[str, Any] | None = None,
        json: Any | None = None,
        method: str = "GET",
    ) -> tuple[dict[str, Any] | str | None, str]:
        """Execute a rate-limited request using ``cfg`` settings."""

        effective_limiter = limiter
        rps = getattr(cfg, "rps", None)
        burst = getattr(cfg, "burst", None)
        if effective_limiter is None and rps:
            effective_limiter = get_limiter(self._resolve_limiter_name(cfg), rps, burst)
        if effective_limiter is not None:
            effective_limiter.acquire()

        delay = self._resolve_delay(cfg)
        timeout = (getattr(cfg, "timeout_connect", 10.0), getattr(cfg, "timeout_read", 10.0))
        retries = getattr(cfg, "retries", 0)

        data, error = _do_request(
            session,
            url,
            delay,
            headers=headers,
            params=params,
            json=json,
            method=method,
            retries=retries,
            timeout=timeout,
            retry_cfg=retry_cfg,
        )
        if error:
            logger.debug("request_fail", extra={"stage": "request_fail", "url": url})
        else:
            logger.debug("request_ok", extra={"stage": "request_ok", "url": url})
        return data, error

