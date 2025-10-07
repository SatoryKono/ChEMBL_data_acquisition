"""HTTP client helpers for interacting with the UniProt REST API."""

from __future__ import annotations

import json
import random
import threading
from collections.abc import Callable, Iterable, Iterator
from contextlib import contextmanager
from typing import Any, cast

import requests
from requests import Session

from ..common.rate_limiter import get_limiter, sleep
from ..config.models import ApiCfg, RetryCfg, UniprotCfg
from ..config.runtime import session_with_retry

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

_session_lock = threading.RLock()
_session_condition = threading.Condition(_session_lock)


def _make_session_factory(api: ApiCfg, retry: RetryCfg) -> Callable[[], Session]:
    def factory() -> Session:
        return session_with_retry(api, retry)

    return factory


_session_factory: Callable[[], Session] = _make_session_factory(ApiCfg(), RetryCfg())
_session_local = threading.local()
_sessions: set[Session] = set()
_active_requests = 0
_retry_cfg: RetryCfg = RetryCfg()
_jitter_provider: Callable[[float], float] | None = None


def _close_sessions(sessions: Iterable[Session]) -> None:
    for session in sessions:
        close = getattr(session, "close", None)
        if callable(close):
            close()


def _get_or_create_session_locked() -> Session:
    session_attr = getattr(_session_local, "session", None)
    session = session_attr if isinstance(session_attr, Session) else None
    if session is None:
        session = _session_factory()
        _sessions.add(session)
        _session_local.session = session
    return session


@contextmanager
def _session_context() -> Iterator[Session]:
    global _active_requests
    with _session_condition:
        session = _get_or_create_session_locked()
        _active_requests += 1
    try:
        yield session
    finally:
        with _session_condition:
            _active_requests -= 1
            if _active_requests == 0:
                _session_condition.notify_all()


def init_session(
    api: ApiCfg,
    retry: RetryCfg,
    *,
    jitter: Callable[[float], float] | None = None,
) -> None:
    """Initialise the shared HTTP session."""

    global _session_factory, _session_local, _retry_cfg, _jitter_provider

    factory = _make_session_factory(api, retry)
    with _session_condition:
        while _active_requests > 0:
            _session_condition.wait()
        old_sessions = list(_sessions)
        _sessions.clear()
        _session_factory = factory
        _session_local = threading.local()
        _retry_cfg = retry
        _jitter_provider = jitter if jitter is not None else retry.build_jitter()

    _close_sessions(old_sessions)


def get_session() -> Session:
    """Expose the configured :class:`requests.Session` instance."""

    with _session_condition:
        return _get_or_create_session_locked()


def fetch_uniprot(uniprot_id: str, *, cfg: UniprotCfg) -> dict[str, Any]:
    """Fetch a UniProt JSON record from the public REST API."""

    base = cfg.base.rstrip("/")
    url = f"{base}/uniprotkb/{uniprot_id}.json"
    timeout = (cfg.timeout_connect, cfg.timeout_read)

    attempt = 1
    while True:
        with _session_condition:
            retry_cfg = _retry_cfg.model_copy(deep=True)
            jitter_provider = _jitter_provider

        if attempt > retry_cfg.max_attempts:
            break

        limiter = get_limiter("uniprot", cfg.rps, cfg.burst)
        limiter.acquire()
        try:
            with _session_context() as session:
                with session.get(url, timeout=timeout) as resp:
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
            if retry_cfg.backoff_factor > 0:
                if jitter_provider is not None:
                    jitter_value = jitter_provider(retry_cfg.backoff_factor)
                else:
                    jitter_value = random.uniform(0, retry_cfg.backoff_factor)
            else:
                jitter_value = 0.0
            delay = backoff + jitter_value + (cfg.delay if cfg.delay else 0)
            if delay > 0:
                sleep(delay)

        attempt += 1

    raise UniProtFetchError(f"UniProt request failed for {uniprot_id}")
