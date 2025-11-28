"""Crossref API helper with resilient retry behaviour."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any, ClassVar

import requests

from library.clients.base import BaseJsonClient
from library.utils.retry import DEFAULT_MAX_TRIES, DEFAULT_TIMEOUT


@dataclass(slots=True)
class CrossrefClient(BaseJsonClient):
    """HTTP client fetching metadata for a given DOI."""

    base_url: str = "https://api.crossref.org/works"
    timeout: float = DEFAULT_TIMEOUT
    max_tries: int = DEFAULT_MAX_TRIES
    session: requests.Session = field(default_factory=requests.Session)
    logger_name: ClassVar[str] = __name__
    client_name: ClassVar[str] = "crossref"
    log_event: ClassVar[str] = "crossref_request"

    def get_metadata(
        self,
        doi: str,
        *,
        params: Mapping[str, Any] | None = None,
        timeout: float | None = None,
    ) -> dict[str, Any]:
        """Return Crossref metadata for ``doi``."""

        endpoint = doi
        effective_timeout = timeout or self.timeout
        return self.request_json(endpoint=endpoint, params=params, timeout=effective_timeout)


__all__ = ["CrossrefClient"]
