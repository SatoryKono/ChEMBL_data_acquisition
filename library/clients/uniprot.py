"""HTTP client helpers for interacting with the UniProt REST API."""

from __future__ import annotations

import json
import random

import threading

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


# Default session uses the application-wide API configuration. Call
# :func:`init_session` with a custom configuration to override it.

_session_lock = threading.Lock()


_session: Session = session_with_retry(ApiCfg(), RetryCfg())
_retry_cfg: RetryCfg = RetryCfg()


def init_session(api: ApiCfg, retry: RetryCfg) -> None:
    """Initialise the shared HTTP session."""

    global _session, _retry_cfg

    new_session = session_with_retry(api, retry)
    old_session: Session | None = None
    with _session_lock:
        old_session = _session
        _session = new_session
        _retry_cfg = retry

    if old_session is not None:
        old_session.close()


def get_session() -> Session:
    """Expose the configured :class:`requests.Session` instance."""

    with _session_lock:
        return _session


def fetch_uniprot(uniprot_id: str, *, cfg: UniprotCfg) -> dict[str, Any]:
    """Fetch a UniProt JSON record from the public REST API."""

    base = cfg.base.rstrip("/")
    url = f"{base}/uniprotkb/{uniprot_id}.json"
    timeout = (cfg.timeout_connect, cfg.timeout_read)

    attempt = 1
    while True:
        with _session_lock:
            retry_cfg = _retry_cfg

        if attempt > retry_cfg.max_attempts:
            break

        limiter = get_limiter("uniprot", cfg.rps, cfg.burst)
        limiter.acquire()
        try:
            with _session_lock:

                with _session.get(url, timeout=timeout) as resp:
                    resp.raise_for_status()
                    try:
                        return cast(dict[str, Any], resp.json())
                    except (
                        json.JSONDecodeError
                    ) as exc:  # pragma: no cover - malformed JSON
                        raise UniProtFetchError(
                            f"Failed to decode JSON for UniProt {uniprot_id}: {exc}"
                        ) from exc

        except requests.RequestException as exc:  # pragma: no cover - network
            if attempt >= retry_cfg.max_attempts:
                raise UniProtFetchError(
                    f"UniProt request failed for {uniprot_id}: {exc}"
                ) from exc

            backoff = retry_cfg.backoff_factor * (2 ** (attempt - 1))
            jitter = random.uniform(0, retry_cfg.backoff_factor)
            delay = backoff + jitter + (cfg.delay if cfg.delay else 0)
            if delay > 0:
                sleep(delay)

        attempt += 1

    raise UniProtFetchError(f"UniProt request failed for {uniprot_id}")
