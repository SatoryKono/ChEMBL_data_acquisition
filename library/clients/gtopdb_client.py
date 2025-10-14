"""Minimal client for Guide to Pharmacology (GtoPdb) services."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import requests

from library.utils.logging import Logger, get_logger
from library.utils.retry import DEFAULT_MAX_TRIES, DEFAULT_TIMEOUT, with_retry


def _build_url(base_url: str, endpoint: str) -> str:
    return f"{base_url.rstrip('/')}/{endpoint.lstrip('/')}"


@dataclass(slots=True)
class GtoPdbClient:
    """HTTP client retrieving GtoPdb target annotations."""

    base_url: str = "https://www.guidetopharmacology.org/services"
    timeout: float = DEFAULT_TIMEOUT
    max_tries: int = DEFAULT_MAX_TRIES
    session: requests.Session = field(default_factory=requests.Session)
    _logger: Logger = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self._logger = get_logger(__name__).bind(client="gtopdb", base_url=self.base_url)

    def get_target_by_uniprot(
        self,
        accession: str,
        *,
        params: Mapping[str, Any] | None = None,
        timeout: float | None = None,
    ) -> Any:
        """Fetch target annotations mapped to ``accession``."""

        endpoint = f"targets/uniprot/{accession}"
        effective_timeout = timeout or self.timeout
        request = with_retry(
            max_tries=self.max_tries,
            timeout=effective_timeout,
            logger=self._logger,
            log_event="gtopdb_request",
        )(self._request_json)
        return request(endpoint=endpoint, params=params, timeout=effective_timeout)

    def _request_json(
        self,
        endpoint: str,
        *,
        params: Mapping[str, Any] | None,
        timeout: float,
    ) -> Any:
        url = _build_url(self.base_url, endpoint)
        self._logger.info("gtopdb_request_start", endpoint=endpoint, url=url)
        response = self.session.get(url, params=params, timeout=timeout)
        response.raise_for_status()
        payload = response.json()
        self._logger.info("gtopdb_request_success", endpoint=endpoint, url=url)
        return payload


__all__ = ["GtoPdbClient"]
