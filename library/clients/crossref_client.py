"""Minimal Crossref client with structured retries."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import requests

from library.utils.logging import StructuredLogger, get_logger, log_context
from library.utils.retry import retryable

__all__ = ["CrossrefServiceClient"]


@dataclass(slots=True)
class CrossrefServiceClient:
    """Access Crossref metadata with resilient HTTP calls."""

    base_url: str
    session: requests.Session = field(default_factory=requests.Session)
    timeout: float = 10.0
    max_attempts: int = 3
    run_id: str | None = None
    logger: StructuredLogger = field(default_factory=lambda: get_logger(__name__))

    def fetch_doi(
        self, doi: str, *, params: Mapping[str, Any] | None = None
    ) -> Any:
        """Return metadata for ``doi``."""

        path = f"works/{doi}"
        return self._get_json(path, params=params)

    def _get_json(
        self, path: str, *, params: Mapping[str, Any] | None = None
    ) -> Any:
        url = f"{self.base_url.rstrip('/')}/{path.lstrip('/')}"

        @retryable(
            logger=self.logger,
            max_attempts=self.max_attempts,
            timeout=self.timeout,
            stage="crossref_request",
        )
        def _send(*, timeout: float) -> requests.Response:
            response = self.session.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response

        with log_context(run_id=self.run_id, stage="crossref_fetch"):
            response = _send()
            payload = response.json()
            self.logger.info(
                "http_success", url=url, status=response.status_code, service="crossref"
            )
            return payload
