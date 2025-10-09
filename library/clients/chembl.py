"""Shared HTTP utilities for ChEMBL API access."""

from __future__ import annotations

import random
import threading
import socket
from collections.abc import Iterable, Iterator
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from itertools import islice
from dataclasses import dataclass, field
from time import monotonic
from types import TracebackType
from typing import Any, Callable, TypeVar, cast
from urllib.parse import urlsplit, urlunsplit

import requests
from cachetools import TTLCache
from requests import Session

try:  # pragma: no cover - urllib3 is always available with requests
    from urllib3.exceptions import ReadTimeoutError as _Urllib3ReadTimeoutError
except Exception:  # pragma: no cover - defensive fallback
    _Urllib3ReadTimeoutError = None  # type: ignore[assignment]

try:  # pragma: no cover - urllib3 is always available with requests
    from urllib3.exceptions import NameResolutionError as _Urllib3NameResolutionError
except Exception:  # pragma: no cover - defensive fallback
    _Urllib3NameResolutionError = None  # type: ignore[assignment]

from ..config.models import ApiCfg, ChemblCacheCfg, RetryCfg
from ..config.runtime import session_with_retry
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
    jitter:
        Optional callable producing jitter values for retry backoff. When not
        provided the jitter is derived from ``retry`` using
        :meth:`library.config.RetryCfg.build_jitter`.
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
        jitter: Callable[[float], float] | None = None,
    ) -> None:
        api = api or ApiCfg()
        retry = retry or RetryCfg()
        self._jitter = jitter if jitter is not None else retry.build_jitter()
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
                    "cache_hit",
                    extra={
                        "url": url,
                        "rps": cfg.rps,
                        "status": "hit",
                        "timeout": read_timeout,
                    },
                )
                return cast(dict[str, Any], cached)
            logger.debug(
                "cache_miss",
                extra={
                    "url": url,
                    "rps": cfg.rps,
                    "status": "miss",
                    "timeout": read_timeout,
                },
            )

        last_exc: BaseException | None = None
        last_exc_cause: BaseException | None = None
        total_attempts = cfg.retries + 1
        fallback_url = _strip_json_suffix(url)

        for attempt in range(1, total_attempts + 1):
            request_url = url
            used_fallback = False
            abort_attempts = False
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
                    extra={
                        "url": request_url,
                        "attempt": attempt,
                        "rps": cfg.rps,
                        "timeout": read_timeout,
                    },
                )
                try:
                    start_time = monotonic()
                    session = self._get_session()
                    with session.get(
                        request_url, timeout=(cfg.timeout_connect, read_timeout)
                    ) as response:
                        response.raise_for_status()
                        try:
                            data = cast(dict[str, Any], response.json())
                        except ValueError as exc:
                            elapsed = monotonic() - start_time
                            logger.exception(
                                "json_error",
                                extra={"url": request_url, "elapsed": elapsed},
                            )
                            raise ValueError(
                                f"invalid JSON in response from {request_url}"
                            ) from exc
                        response_elapsed = getattr(response, "elapsed", None)
                        response_elapsed_seconds: float | None
                        if (
                            response_elapsed is not None
                            and hasattr(response_elapsed, "total_seconds")
                        ):
                            response_elapsed_seconds = float(
                                response_elapsed.total_seconds()
                            )
                        else:
                            response_elapsed_seconds = None
                        duration = monotonic() - start_time
                        logger.debug(
                            "request_ok",
                            extra={
                                "url": request_url,
                                "status": getattr(response, "status_code", None),
                                "rps": cfg.rps,
                                "elapsed": duration,
                                "response_elapsed": response_elapsed_seconds,
                                "timeout": read_timeout,
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
                                    "elapsed": duration,
                                    "response_elapsed": response_elapsed_seconds,
                                    "timeout": read_timeout,
                                },
                            )
                        return data
                except ValueError as exc:
                    elapsed = monotonic() - start_time
                    last_exc = exc
                    last_exc_cause = None
                    if attempt >= total_attempts:
                        logger.exception(
                            "request_fail",
                            extra={
                                "url": request_url,
                                "status": None,
                                "rps": cfg.rps,
                                "elapsed": elapsed,
                                "attempt": attempt,
                                "timeout": read_timeout,
                                "retry": attempt,
                                "backoff": 0.0,
                            },
                        )
                        break
                    delay = _backoff_delay(attempt, cfg, header_delay=None, jitter=self._jitter)
                    _log_retry_warning(
                        "request_retry_json_error",
                        url=request_url,
                        attempt=attempt,
                        delay=delay,
                        elapsed=elapsed,
                        status=None,
                        exc_info=True,
                    )
                    _log_retry_delay(request_url, attempt, None, delay)
                    sleep(delay)
                    break
                except requests.HTTPError as exc:
                    elapsed = monotonic() - start_time
                    last_exc = exc
                    last_exc_cause = None
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
                                "elapsed": elapsed,
                            },
                        )
                        continue
                    if attempt >= total_attempts:
                        logger.exception(
                            "request_fail",
                            extra={
                                "url": request_url,
                                "status": status,
                                "rps": cfg.rps,
                                "elapsed": elapsed,
                                "attempt": attempt,
                                "timeout": read_timeout,
                                "retry": attempt,
                                "backoff": 0.0,
                            },
                        )
                        break
                    header_delay = _retry_after_delay(response)
                    delay = _backoff_delay(
                        attempt, cfg, header_delay, jitter=self._jitter
                    )
                    _log_retry_warning(
                        "request_retry_http_error",
                        url=request_url,
                        attempt=attempt,
                        delay=delay,
                        elapsed=elapsed,
                        status=status,
                        exc_info=True,
                    )
                    _log_retry_delay(request_url, attempt, status, delay, header_delay)
                    sleep(delay)
                    break
                except requests.RequestException as exc:
                    elapsed = monotonic() - start_time
                    normalized_exc = _normalise_request_exception(exc)
                    last_exc = normalized_exc
                    last_exc_cause = exc if normalized_exc is not exc else None
                    response = getattr(exc, "response", None)
                    status = getattr(response, "status_code", None)
                    name_resolution_error = _is_name_resolution_error(exc)
                    if (
                        not used_fallback
                        and fallback_url is not None
                        and fallback_url != request_url
                        and _should_switch_to_fallback(normalized_exc)
                        and not name_resolution_error
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
                                "elapsed": elapsed,
                            },
                        )
                        continue
                    if name_resolution_error:
                        logger.exception(
                            "request_name_resolution_error",
                            extra={
                                "url": request_url,
                                "attempt": attempt,
                                "rps": cfg.rps,
                                "timeout": read_timeout,
                                "retry": attempt,
                                "backoff": 0.0,
                                "status": status,
                                "elapsed": elapsed,
                            },
                        )
                        abort_attempts = True
                        break
                    if attempt >= total_attempts:
                        logger.exception(
                            "request_fail",
                            extra={
                                "url": request_url,
                                "status": None,
                                "rps": cfg.rps,
                                "elapsed": elapsed,
                                "attempt": attempt,
                                "timeout": read_timeout,
                                "retry": attempt,
                                "backoff": 0.0,
                            },
                        )
                        break
                    delay = _backoff_delay(
                        attempt, cfg, header_delay=None, jitter=self._jitter
                    )
                    _log_retry_warning(
                        "request_retry_exception",
                        url=request_url,
                        attempt=attempt,
                        delay=delay,
                        elapsed=elapsed,
                        status=status,
                        exc_info=True,
                    )
                    _log_retry_delay(request_url, attempt, status, delay)
                    sleep(delay)
                    break

            if abort_attempts:
                break

        if last_exc is not None:
            if last_exc_cause is not None and last_exc is not last_exc_cause:
                raise last_exc from last_exc_cause
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
    attempt: int,
    cfg: ApiCfg,
    header_delay: float | None,
    *,
    jitter: Callable[[float], float] | None = None,
) -> float:
    base = cfg.backoff_factor * (2 ** (attempt - 1))
    jitter_value = 0.0
    if cfg.backoff_factor > 0:
        if jitter is not None:
            jitter_value = jitter(base)
        else:
            jitter_value = random.uniform(0, base)
    delay = base + jitter_value
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


def _log_retry_warning(
    event: str,
    *,
    url: str,
    attempt: int,
    delay: float,
    elapsed: float,
    status: int | None,
    exc_info: bool = False,
) -> None:
    logger.warning(
        event,
        extra={
            "url": url,
            "retry": attempt,
            "backoff": delay,
            "status": status,
            "elapsed": elapsed,
        },
        exc_info=exc_info,
    )


def _should_switch_to_fallback(exception: requests.RequestException) -> bool:
    """Return ``True`` when ``exception`` warrants trying the fallback URL."""

    return isinstance(
        exception,
        (
            requests.Timeout,
            requests.ConnectionError,
        ),
    )


def _normalise_request_exception(
    exc: requests.RequestException,
) -> requests.RequestException:
    """Return a timeout-aware variant of ``exc`` when possible."""

    if isinstance(exc, requests.ReadTimeout):
        return exc

    timeout_error = _find_read_timeout_error(exc)
    if timeout_error is not None:
        message = str(timeout_error) or "Read timed out."
        request = getattr(timeout_error, "request", None)
        if request is None:
            request = getattr(exc, "request", None)
        response = getattr(timeout_error, "response", None)
        if response is None:
            response = getattr(exc, "response", None)
        return requests.ReadTimeout(
            message,
            request=request,
            response=response,
        )

    if isinstance(exc, requests.ConnectionError):
        message = str(exc)
        lowered = message.lower()
        if "read" in lowered and ("timeout" in lowered or "timed out" in lowered):
            normalized_message = message or "Read timed out."
            return requests.ReadTimeout(
                normalized_message,
                request=getattr(exc, "request", None),
                response=getattr(exc, "response", None),
            )

    return exc


def _find_read_timeout_error(exc: BaseException) -> BaseException | None:
    """Return the first read-timeout error found in ``exc``'s chain."""

    seen_ids: set[int] = set()
    current: BaseException | None = exc
    timeout_types: tuple[type[BaseException], ...]
    if _Urllib3ReadTimeoutError is None:  # pragma: no cover - defensive
        timeout_types = (requests.ReadTimeout,)
    else:
        timeout_types = (requests.ReadTimeout, _Urllib3ReadTimeoutError)

    while current is not None:
        current_id = id(current)
        if current_id in seen_ids:
            return None
        seen_ids.add(current_id)
        if isinstance(current, timeout_types):
            return current
        for arg in getattr(current, "args", ()):
            if isinstance(arg, BaseException) and id(arg) not in seen_ids:
                nested = _find_read_timeout_error(arg)
                if nested is not None:
                    return nested
        reason = getattr(current, "reason", None)
        if isinstance(reason, BaseException) and id(reason) not in seen_ids:
            nested_reason = _find_read_timeout_error(reason)
            if nested_reason is not None:
                return nested_reason
        cause = current.__cause__
        if isinstance(cause, BaseException):
            current = cause
            continue
        context = current.__context__
        current = context if isinstance(context, BaseException) else None
    return None


def _iter_exception_chain(exc: BaseException) -> Iterator[BaseException]:
    """Yield ``exc`` and nested exceptions including causes and contexts."""

    seen: set[int] = set()
    stack: list[BaseException] = [exc]
    while stack:
        candidate = stack.pop()
        candidate_id = id(candidate)
        if candidate_id in seen:
            continue
        seen.add(candidate_id)
        yield candidate
        cause = candidate.__cause__
        if isinstance(cause, BaseException):
            stack.append(cause)
        context = candidate.__context__
        if isinstance(context, BaseException):
            stack.append(context)
        reason = getattr(candidate, "reason", None)
        if isinstance(reason, BaseException):
            stack.append(reason)
        for argument in getattr(candidate, "args", ()):  # pragma: no branch - tuple walk
            if isinstance(argument, BaseException):
                stack.append(argument)


def _is_name_resolution_error(exc: BaseException) -> bool:
    """Return ``True`` when ``exc`` represents a DNS resolution failure."""

    indicators = ("name resolution", "nameresolutionerror", "getaddrinfo failed")
    for candidate in _iter_exception_chain(exc):
        if (
            _Urllib3NameResolutionError is not None
            and isinstance(candidate, _Urllib3NameResolutionError)
        ):
            return True
        if isinstance(candidate, socket.gaierror):
            return True
        message = str(candidate).strip().lower()
        if any(indicator in message for indicator in indicators):
            return True
    return False


__all__ = ["ChemblClient", "_chunked"]
