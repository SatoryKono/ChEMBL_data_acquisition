"""Retry utilities ensuring resilient HTTP calls."""

from __future__ import annotations

import functools
from typing import Any, Callable, Iterator, ParamSpec, TypeVar

import requests

from .logging import Logger, get_logger

try:  # pragma: no cover - prefer real backoff when available
    import backoff
except ModuleNotFoundError:  # pragma: no cover - fallback for minimal environments
    class _FallbackBackoff:
        @staticmethod
        def expo(base: float = 1.0, factor: float = 2.0, max_value: float | None = None):
            def _generator() -> Iterator[float]:
                delay = base
                while True:
                    yield delay if max_value is None else min(delay, max_value)
                    delay *= factor

            return _generator()

        @staticmethod
        def full_jitter(value: float) -> float:
            return value

        @staticmethod
        def on_exception(
            wait_gen: Callable[..., Iterator[float]],
            exception: type[BaseException] | tuple[type[BaseException], ...],
            *,
            max_tries: int = 1,
            jitter: Callable[[float], float] | None = None,
            on_backoff: Callable[[dict[str, Any]], None] | None = None,
            on_giveup: Callable[[dict[str, Any]], None] | None = None,
        ) -> Callable[[Callable[P, R]], Callable[P, R]]:
            def decorator(func: Callable[P, R]) -> Callable[P, R]:
                def wrapped(*args: P.args, **kwargs: P.kwargs) -> R:
                    tries = 0
                    delays = wait_gen()
                    while True:
                        try:
                            return func(*args, **kwargs)
                        except exception as exc:  # type: ignore[misc]
                            tries += 1
                            details = {
                                "target": func,
                                "args": args,
                                "kwargs": kwargs,
                                "tries": tries,
                                "value": exc,
                            }
                            if tries >= max_tries:
                                if on_giveup is not None:
                                    on_giveup(details)
                                raise
                            wait = next(delays)
                            if jitter is not None:
                                wait = jitter(wait)
                            if on_backoff is not None:
                                details["wait"] = wait
                                on_backoff(details)
                            continue

                return wrapped

            return decorator

    backoff = _FallbackBackoff()

P = ParamSpec("P")
R = TypeVar("R")

DEFAULT_TIMEOUT = 10.0
DEFAULT_MAX_TRIES = 3


def _bind_logger(logger: Logger | None, module_name: str, event: str) -> Logger:
    base = logger if logger is not None else get_logger(module_name)
    return base.bind(event=event)


def with_retry(
    *,
    max_tries: int = DEFAULT_MAX_TRIES,
    timeout: float = DEFAULT_TIMEOUT,
    logger: Logger | None = None,
    log_event: str = "http_request",
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Return a decorator applying exponential backoff to ``func``.

    Parameters
    ----------
    max_tries:
        Maximum number of attempts including the first call.
    timeout:
        Default timeout (seconds) enforced when the wrapped function does not
        specify one explicitly.
    logger:
        Optional structured logger used for retry diagnostics. When omitted a
        module-level logger is created using :func:`get_logger`.
    log_event:
        Base event name used when emitting retry logs.
    """

    if max_tries < 1:
        raise ValueError("max_tries must be at least 1")
    if timeout <= 0:
        raise ValueError("timeout must be positive")

    def decorator(func: Callable[P, R]) -> Callable[P, R]:
        bound_logger = _bind_logger(logger, func.__module__, log_event)

        def _on_backoff(details: dict[str, Any]) -> None:
            attempt = details.get("tries", 0)
            wait = details.get("wait")
            exc = details.get("value")
            kwargs = details.get("kwargs", {})
            payload: dict[str, Any] = {
                "attempt": attempt,
                "max_tries": max_tries,
            }
            endpoint = kwargs.get("endpoint")
            if endpoint is not None:
                payload["endpoint"] = endpoint
            if wait is not None:
                payload["sleep_s"] = wait
            if isinstance(exc, Exception):
                payload["exception"] = exc.__class__.__name__
                payload["error"] = str(exc)
            bound_logger.warning(f"{log_event}_retry", **payload)

        def _on_giveup(details: dict[str, Any]) -> None:
            exc = details.get("value")
            payload: dict[str, Any] = {"attempts": details.get("tries", max_tries)}
            kwargs = details.get("kwargs", {})
            endpoint = kwargs.get("endpoint")
            if endpoint is not None:
                payload["endpoint"] = endpoint
            if isinstance(exc, Exception):
                payload["exception"] = exc.__class__.__name__
                payload["error"] = str(exc)
            bound_logger.error(f"{log_event}_giveup", **payload)

        @backoff.on_exception(
            backoff.expo,
            requests.exceptions.RequestException,
            max_tries=max_tries,
            jitter=backoff.full_jitter,
            on_backoff=_on_backoff,
            on_giveup=_on_giveup,
        )
        @functools.wraps(func)
        def wrapped(*args: P.args, **kwargs: P.kwargs) -> R:
            kwargs.setdefault("timeout", timeout)
            return func(*args, **kwargs)

        return wrapped

    return decorator


__all__ = ["DEFAULT_MAX_TRIES", "DEFAULT_TIMEOUT", "with_retry"]
