"""PubChem API helper with deterministic retries."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any, ClassVar

import requests

from library.clients.base import BaseJsonClient
from library.utils.retry import DEFAULT_MAX_TRIES, DEFAULT_TIMEOUT


@dataclass(slots=True)
class PubChemClient(BaseJsonClient):
    """HTTP client exposing the PubChem compound endpoint."""

    base_url: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    timeout: float = DEFAULT_TIMEOUT
    max_tries: int = DEFAULT_MAX_TRIES
    session: requests.Session = field(default_factory=requests.Session)
    logger_name: ClassVar[str] = __name__
    client_name: ClassVar[str] = "pubchem"
    log_event: ClassVar[str] = "pubchem_request"

    def get_compound(
        self,
        cid: str,
        *,
        params: Mapping[str, Any] | None = None,
        timeout: float | None = None,
    ) -> dict[str, Any]:
        """Fetch compound properties for ``cid``."""

        endpoint = f"compound/cid/{cid}/JSON"
        effective_timeout = timeout or self.timeout
        return self.request_json(endpoint=endpoint, params=params, timeout=effective_timeout)


__all__ = ["PubChemClient"]
