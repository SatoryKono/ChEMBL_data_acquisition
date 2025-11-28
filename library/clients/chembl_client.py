"""Lightweight Chembl API client with structured logging."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, ClassVar

import requests

from library.clients.base import BaseJsonClient
from library.utils.retry import DEFAULT_MAX_TRIES, DEFAULT_TIMEOUT


@dataclass(slots=True)
class ChemblClient(BaseJsonClient):
    """HTTP client providing minimal Chembl read access."""

    base_url: str = "https://www.ebi.ac.uk/chembl/api/data"
    timeout: float = DEFAULT_TIMEOUT
    max_tries: int = DEFAULT_MAX_TRIES
    session: requests.Session = field(default_factory=requests.Session)
    logger_name: ClassVar[str] = __name__
    client_name: ClassVar[str] = "chembl"
    log_event: ClassVar[str] = "chembl_request"

    def get_molecule(
        self,
        chembl_id: str,
        *,
        params: Mapping[str, Any] | None = None,
        timeout: float | None = None,
    ) -> dict[str, Any]:
        """Fetch molecule metadata for ``chembl_id``."""

        endpoint = f"molecule/{chembl_id}"
        effective_timeout = timeout or self.timeout
        return self.request_json(endpoint=endpoint, params=params, timeout=effective_timeout)

    def list_targets(
        self,
        *,
        limit: int,
        offset: int = 0,
        fields: Sequence[str] | None = None,
        params: Mapping[str, Any] | None = None,
        timeout: float | None = None,
    ) -> dict[str, Any]:
        """Return a page of targets with optional field selection."""

        if limit <= 0:
            raise ValueError("limit must be positive")

        effective_params: dict[str, Any] = {"limit": limit, "offset": offset}
        if fields:
            effective_params["fields"] = ",".join(fields)
        if params:
            effective_params.update(params)

        endpoint = "target.json"
        effective_timeout = timeout or self.timeout
        return self.request_json(
            endpoint=endpoint, params=effective_params, timeout=effective_timeout
        )


__all__ = ["ChemblClient"]
