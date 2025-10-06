"""Shared HTTP utilities for ChEMBL API access."""

from __future__ import annotations

import random
import threading
from collections.abc import Iterable, Iterator
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from itertools import islice
from dataclasses import dataclass, field
from time import perf_counter
from types import TracebackType
from typing import Any, Callable, TypeVar, cast
from urllib.parse import urlsplit, urlunsplit

import requests
from cachetools import TTLCache
from requests import Session

from ..config import ApiCfg, ChemblCacheCfg, RetryCfg, session_with_retry
from ..common.log import logger
from ..common.rate_limiter import RateLimiter, get_limiter, sleep


@dataclass
class ChemblClient:
    """HTTP client for the ChEMBL API with a TTL cache.

    Parameters
    ----------
    api:
        Global API settings providing the ``User-Agent`` header.
    retry:
        Retry configuration applied to all requests.
    chembl:
        Optional ChEMBL-specific configuration controlling cache TTL and size.
    session:
        Optional pre-configured :class:`requests.Session` instance; primarily
        intended for tests.
    global_limiter:
        Optional system-wide :class:`RateLimiter` enforcing ``Config.rate``
        across all HTTP clients.
    """

    cache: TTLCache[str, dict[str, Any]] = field(init=False)
    _cache_lock: threading.Lock = field(default_factory=threading.Lock, init=False)
    _session_local: threading.local = field(init=False)
    _sessions: set[Session] = field(init=False)
    _sessions_lock: threading.Lock = field(default_factory=threading.Lock, init=False)
    _session_factory: Callable[[], Session] = field(init=False)
    _global_limiter: RateLimiter | None = field(default=None, init=False)

    def __init__(
        self,
        api: ApiCfg | None = None,
        retry: RetryCfg | None = None,
        chembl: ChemblCacheCfg | None = None,
        *,
        session: Session | None = None,
        global_limiter: RateLimiter | None = None,
    ) -> None:
        api = api or ApiCfg()
        retry = retry or RetryCfg()
        if session is not None:
            def _session_from_argument(provided: Session = session) -> Session:
                return provided

            self._session_factory = _session_from_argument
        else:
            api_cfg_default: ApiCfg = api
            retry_cfg_default: RetryCfg = retry

            def _build_session(
                api_cfg: ApiCfg = api_cfg_default,
                retry_cfg: RetryCfg = retry_cfg_default,
            ) -> Session:
                return session_with_retry(api_cfg, retry_cfg)

            self._session_factory = _build_session
        ttl = chembl.cache_ttl if chembl is not None else ChemblCacheCfg().cache_ttl
        maxsize = (
            chembl.cache_maxsize
            if chembl is not None
            else ChemblCacheCfg().cache_maxsize
        )
        self.cache = TTLCache(maxsize=maxsize, ttl=ttl)
        self._cache_lock = threading.Lock()
        self._session_local = threading.local()
        self._sessions = set()
        self._sessions_lock = threading.Lock()
        self._global_limiter = global_limiter

    def close(self) -> None:
        """Close the underlying HTTP session.

        This method releases network resources held by the internal
        :class:`requests.Session` instance. It is safe to call multiple times.
        """

        with self._sessions_lock:
            sessions = list(self._sessions)
            self._sessions.clear()
        for session in sessions:
            close = getattr(session, "close", None)
            if callable(close):
                close()
        self._session_local = threading.local()

    def __enter__(self) -> ChemblClient:
        """Return ``self`` when entering a context manager."""

        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        """Close the session when leaving a context manager.

        Parameters
        ----------
        exc_type, exc, tb:
            Exception information supplied by the runtime; unused.
        """

        self.close()
        return None

    def _register_session(self, session: Session) -> None:
        with self._sessions_lock:
            self._sessions.add(session)

    def _get_session(self) -> Session:
        session_attr = getattr(self._session_local, "session", None)
        session = (
            session_attr
            if isinstance(session_attr, Session)
            or (hasattr(session_attr, "get") and hasattr(session_attr, "close"))
            else None
        )
        if session is None:
            session = self._session_factory()
            self._register_session(session)
            setattr(self._session_local, "session", session)
        return cast(Session, session)

    @property
    def session(self) -> Session:
        """Return the session bound to the current thread."""

        return self._get_session()

    def request_json(
        self, url: str, *, cfg: ApiCfg, timeout: float | None = None
    ) -> dict[str, Any]:
        """Return JSON content from ``url``.

        Parameters
        ----------
        url:
            API endpoint to query.
        cfg:
            Configuration providing timeout and retry settings.
        timeout:
            Optional override for the read timeout in seconds.

        Notes
        -----
        The response body is decoded using the declared encoding or UTF-8,
        replacing undecodable bytes with ``\ufffd`` before JSON parsing.

        Returns
        -------
        dict[str, Any]
            Parsed JSON document.

        Raises
        ------
        requests.RequestException
            If the HTTP request fails.
        ValueError
            If the response body is not valid JSON or cannot be decoded.
        """

        limiter = get_limiter("chembl", cfg.rps, cfg.burst)
        read_timeout = timeout if timeout is not None else cfg.timeout_read
        cache_key = url
        with self._cache_lock:
            cached = self.cache.get(cache_key)
            if cached is not None:
                logger.debug(
                    "cache_hit", extra={"url": url, "rps": cfg.rps, "status": "hit"}
                )
                return cast(dict[str, Any], cached)
            logger.debug(
                "cache_miss", extra={"url": url, "rps": cfg.rps, "status": "miss"}
            )

        last_exc: requests.RequestException | ValueError | None = None
        total_attempts = cfg.retries + 1
        fallback_url = _strip_json_suffix(url)

        for attempt in range(1, total_attempts + 1):
            request_url = url
            used_fallback = False
            while True:
                if self._global_limiter is not None:
                    self._global_limiter.acquire()
                limiter.acquire()
                if used_fallback:
                    event = "request_fallback"
                else:
                    event = "request_start" if attempt == 1 else "request_retry"
                logger.debug(
                    event,
                    extra={"url": request_url, "attempt": attempt, "rps": cfg.rps},
                )
                try:
                    session = self._get_session()
                    start_time = perf_counter()
                    with session.get(
                        request_url, timeout=(cfg.timeout_connect, read_timeout)
                    ) as response:
                        response.raise_for_status()
                        try:
                            data = cast(dict[str, Any], response.json())
                        except ValueError as exc:
                            logger.exception("json_error", extra={"url": request_url})
                            raise ValueError(
                                f"invalid JSON in response from {request_url}"
                            ) from exc
                        elapsed = getattr(response, "elapsed", None)
                        if elapsed is not None and hasattr(elapsed, "total_seconds"):
                            duration = elapsed.total_seconds()
                        else:
                            duration = perf_counter() - start_time
                        logger.debug(
                            "request_ok",
                            extra={
                                "url": request_url,
                                "status": getattr(response, "status_code", None),
                                "rps": cfg.rps,
                                "elapsed": duration,
                            },
                        )
                        with self._cache_lock:
                            cached = self.cache.get(cache_key)
                            if cached is not None:
                                return cast(dict[str, Any], cached)
                            self.cache[cache_key] = data
                            logger.debug(
                                "cache_set", extra={"url": url, "rps": cfg.rps}
                            )
                        if used_fallback:
                            logger.debug(
                                "request_fallback_ok",
                                extra={
                                    "original_url": url,
                                    "fallback_url": request_url,
                                    "attempt": attempt,
                                    "rps": cfg.rps,
                                },
                            )
                        return data
                except ValueError as exc:
                    last_exc = exc
                    if attempt >= total_attempts:
                        logger.exception(
                            "request_fail",
                            extra={"url": request_url, "status": None, "rps": cfg.rps},
                        )
                        break
                    delay = _backoff_delay(attempt, cfg, header_delay=None)
                    _log_retry_delay(request_url, attempt, None, delay)
                    sleep(delay)
                    break
                except requests.HTTPError as exc:
                    last_exc = exc
                    response = exc.response
                    status = getattr(response, "status_code", None)
                    if (
                        status == 404
                        and not used_fallback
                        and fallback_url is not None
                        and fallback_url != request_url
                    ):
                        used_fallback = True
                        request_url = fallback_url
                        logger.debug(
                            "request_fallback_switch",
                            extra={
                                "original_url": url,
                                "fallback_url": request_url,
                                "attempt": attempt,
                                "rps": cfg.rps,
                            },
                        )
                        continue
                    if attempt >= total_attempts:
                        logger.exception(
                            "request_fail",
                            extra={"url": request_url, "status": status, "rps": cfg.rps},
                        )
                        break
                    header_delay = _retry_after_delay(response)
                    delay = _backoff_delay(attempt, cfg, header_delay)
                    _log_retry_delay(request_url, attempt, status, delay, header_delay)
                    sleep(delay)
                    break
                except requests.RequestException as exc:
                    last_exc = exc
                    if attempt >= total_attempts:
                        logger.exception(
                            "request_fail",
                            extra={"url": request_url, "status": None, "rps": cfg.rps},
                        )
                        break
                    delay = _backoff_delay(attempt, cfg, header_delay=None)
                    _log_retry_delay(request_url, attempt, None, delay)
                    sleep(delay)
                    break

        if last_exc is not None:
            raise last_exc
        raise RuntimeError(f"Request loop exited unexpectedly for {url}")

    def clear_cache(self) -> None:
        """Remove all entries from the in-memory cache."""

        with self._cache_lock:
            self.cache.clear()


T = TypeVar("T")


def _chunked(items: Iterable[T], size: int) -> Iterator[list[T]]:
    """Yield ``size``-sized lists from *items*.

    Parameters
    ----------
    items:
        Iterable of identifiers to split. Accepts generators and other lazy
        iterables.
    size:
        Desired chunk size; must be positive.

    Yields
    ------
    list[str]
        Subsequences of ``items`` with at most ``size`` elements.

    Raises
    ------
    ValueError
        If ``size`` is not a positive integer.
    """

    if size <= 0:
        raise ValueError("size must be a positive integer")

    iterator = iter(items)
    while True:
        chunk = list(islice(iterator, size))
        if not chunk:
            break
        yield chunk


def _strip_json_suffix(url: str) -> str | None:
    """Return ``url`` without a trailing ``.json`` path component if present."""

    split = urlsplit(url)
    path = split.path
    if not path.endswith(".json"):
        return None
    stripped = path[:-5]
    if not stripped:
        return None
    return urlunsplit(split._replace(path=stripped))


def _utcnow() -> datetime:
    return datetime.now(timezone.utc)


def _retry_after_delay(response: requests.Response | None) -> float | None:
    if response is None:
        return None
    status = getattr(response, "status_code", None)
    if status is None or not _is_retry_after_applicable(status):
        return None
    header = response.headers.get("Retry-After") if hasattr(response, "headers") else None
    if header is None:
        return None
    value = header.strip()
    if not value:
        return None
    try:
        delay = float(value)
        return max(0.0, delay)
    except ValueError:
        try:
            parsed_dt = parsedate_to_datetime(value)
        except (TypeError, ValueError):
            return None
        if parsed_dt is None:
            return None
        dt = parsed_dt
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        delta = float((dt - _utcnow()).total_seconds())
        return max(0.0, delta)


def _is_retry_after_applicable(status: int) -> bool:
    return status == 429 or 500 <= status < 600


def _backoff_delay(
    attempt: int, cfg: ApiCfg, header_delay: float | None
) -> float:
    base = cfg.backoff_factor * (2 ** (attempt - 1))
    jitter = random.uniform(0, cfg.backoff_factor)
    delay = base + jitter
    if header_delay is not None:
        delay = max(delay, header_delay)
    return float(delay)


def _log_retry_delay(
    url: str,
    attempt: int,
    status: int | None,
    delay: float,
    header_delay: float | None = None,
) -> None:
    logger.debug(
        "retry_sleep",
        extra={
            "url": url,
            "attempt": attempt,
            "status": status,
            "delay": delay,
            "retry_after": header_delay,
        },
    )


__all__ = ["ChemblClient", "_chunked"]
