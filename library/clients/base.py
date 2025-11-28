"""Shared helpers for JSON-based HTTP clients."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any, ClassVar

import requests

from library.utils.logging import Logger, get_logger
from library.utils.retry import DEFAULT_MAX_TRIES, DEFAULT_TIMEOUT, with_retry


def build_url(base_url: str, endpoint: str) -> str:
    """Return ``endpoint`` joined to ``base_url`` with a single ``/`` separator."""

    return f"{base_url.rstrip('/')}/{endpoint.lstrip('/')}"


@dataclass(slots=True)
class BaseJsonClient:
    """Minimal reusable JSON HTTP client with retry logging."""

    base_url: str = ""
    timeout: float = DEFAULT_TIMEOUT
    max_tries: int = DEFAULT_MAX_TRIES
    session: requests.Session = field(default_factory=requests.Session)
    logger: Logger = field(init=False, repr=False)

    logger_name: ClassVar[str] = __name__
    client_name: ClassVar[str] = "base_json_client"
    log_event: ClassVar[str] = "json_request"

    def __post_init__(self) -> None:
        self.logger = get_logger(self.logger_name).bind(
            client=self.client_name, base_url=self.base_url
        )

    def request_json(
        self,
        endpoint: str,
        *,
        params: Mapping[str, Any] | None = None,
        timeout: float | None = None,
    ) -> Any:
        """Issue a GET request and parse the JSON payload with retries."""

        effective_timeout = timeout or self.timeout
        request = with_retry(
            max_tries=self.max_tries,
            timeout=effective_timeout,
            logger=self.logger,
            log_event=self.log_event,
        )(self._request_json)
        return request(endpoint=endpoint, params=params, timeout=effective_timeout)

    def _request_json(
        self,
        endpoint: str,
        *,
        params: Mapping[str, Any] | None,
        timeout: float,
    ) -> Any:
        url = build_url(self.base_url, endpoint)
        self.logger.info(f"{self.log_event}_start", endpoint=endpoint, url=url)
        response = self.session.get(url, params=params, timeout=timeout)
        response.raise_for_status()
        payload = response.json()
        self.logger.info(f"{self.log_event}_success", endpoint=endpoint, url=url)
        return payload


__all__ = ["BaseJsonClient", "build_url"]
