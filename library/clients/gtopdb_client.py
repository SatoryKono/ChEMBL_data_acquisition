"""Minimal client for Guide to Pharmacology (GtoPdb) services."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, ClassVar, Mapping

import requests

from library.clients.base import BaseJsonClient
from library.utils.retry import DEFAULT_MAX_TRIES, DEFAULT_TIMEOUT


@dataclass(slots=True)
class GtoPdbClient(BaseJsonClient):
    """HTTP client retrieving GtoPdb target annotations."""

    base_url: str = "https://www.guidetopharmacology.org/services"
    timeout: float = DEFAULT_TIMEOUT
    max_tries: int = DEFAULT_MAX_TRIES
    session: requests.Session = field(default_factory=requests.Session)
    logger_name: ClassVar[str] = __name__
    client_name: ClassVar[str] = "gtopdb"
    log_event: ClassVar[str] = "gtopdb_request"

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
        return self.request_json(endpoint=endpoint, params=params, timeout=effective_timeout)


__all__ = ["GtoPdbClient"]
