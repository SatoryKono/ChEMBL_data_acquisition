"""HTTP client helpers for interacting with the UniProt REST API."""

from __future__ import annotations

import json
import random
from typing import Any, cast

import requests
from requests import Session

from ..config import ApiCfg, RetryCfg, UniprotCfg, session_with_retry
from ..rate_limiter import get_limiter, sleep

__all__ = [
    "UniProtFetchError",
    "init_session",
    "fetch_uniprot",
    "get_session",
]


class UniProtFetchError(RuntimeError):
    """Raised when a UniProt record cannot be retrieved or decoded."""


# Default session using placeholder contact details. Call :func:`init_session`
# with a proper configuration to set your own user agent.
_session: Session = session_with_retry(
    ApiCfg(user_agent="chembl-da/0.1 (mailto:contact@example.org)"), RetryCfg()
)
_retry_cfg: RetryCfg = RetryCfg()


def init_session(api: ApiCfg, retry: RetryCfg) -> None:
    """Initialise the shared HTTP session."""

    global _session, _retry_cfg
    _session = session_with_retry(api, retry)
    _retry_cfg = retry


def get_session() -> Session:
    """Expose the configured :class:`requests.Session` instance."""

    return _session


def fetch_uniprot(uniprot_id: str, *, cfg: UniprotCfg) -> dict[str, Any]:
    """Fetch a UniProt JSON record from the public REST API."""

    base = cfg.base.rstrip("/")
    url = f"{base}/uniprotkb/{uniprot_id}.json"
    timeout = (cfg.timeout_connect, cfg.timeout_read)

    for attempt in range(1, _retry_cfg.max_attempts + 1):
        limiter = get_limiter("uniprot", cfg.rps, cfg.burst)
        limiter.acquire()
        try:
            with _session.get(url, timeout=timeout) as resp:
                resp.raise_for_status()
                try:
                    return cast(dict[str, Any], resp.json())
                except json.JSONDecodeError as exc:  # pragma: no cover - malformed JSON
                    raise UniProtFetchError(
                        f"Failed to decode JSON for UniProt {uniprot_id}: {exc}"
                    ) from exc
        except requests.RequestException as exc:  # pragma: no cover - network
            if attempt >= _retry_cfg.max_attempts:
                raise UniProtFetchError(
                    f"UniProt request failed for {uniprot_id}: {exc}"
                ) from exc

            backoff = _retry_cfg.backoff_factor * (2 ** (attempt - 1))
            jitter = random.uniform(0, _retry_cfg.backoff_factor)
            delay = backoff + jitter + (cfg.delay if cfg.delay else 0)
            if delay > 0:
                sleep(delay)

    raise UniProtFetchError(f"UniProt request failed for {uniprot_id}")
