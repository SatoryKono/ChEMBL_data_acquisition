"""Reusable retry helpers for resilient HTTP calls."""

from __future__ import annotations

import functools
import time
from typing import Any, Callable, Iterator, TypeVar

try:  # pragma: no cover - production installs provide backoff
    import backoff  # type: ignore[import-not-found]
except ModuleNotFoundError:  # pragma: no cover - fallback for minimal environments
    class _FallbackBackoff:
        @staticmethod
        def expo(**_: Any) -> Callable[[], Iterator[float]]:
            def _generator() -> Iterator[float]:
                delay = 0.0
                while True:
                    yield delay
                    delay = delay * 2 if delay else 1.0

            return _generator()

        @staticmethod
        def on_exception(
            wait_gen: Callable[[], Iterator[float]],
            exception: type[BaseException] | tuple[type[BaseException], ...],
            *,
            max_tries: int | None = None,
            on_backoff: Callable[[dict[str, Any]], None] | None = None,
            on_giveup: Callable[[dict[str, Any]], None] | None = None,
        ) -> Callable[[F], F]:
            def decorator(func: F) -> F:
                @functools.wraps(func)
                def wrapper(*args: Any, **kwargs: Any) -> Any:
                    tries = 0
                    waits = wait_gen()
                    while True:
                        try:
                            return func(*args, **kwargs)
                        except exception as exc:  # type: ignore[misc]
                            tries += 1
                            wait = next(waits, 0.0)
                            payload = {
                                "exception": exc,
                                "tries": tries,
                                "wait": wait,
                                "args": args,
                                "kwargs": kwargs,
                            }
                            if max_tries is not None and tries >= max_tries:
                                if on_giveup:
                                    on_giveup(payload)
                                raise
                            if on_backoff:
                                on_backoff(payload)
                            if wait > 0:
                                time.sleep(wait)

                return wrapper  # type: ignore[return-value]

            return decorator

    backoff = _FallbackBackoff()

import requests

from .logging import StructuredLogger, get_logger

__all__ = ["DEFAULT_MAX_ATTEMPTS", "DEFAULT_TIMEOUT", "retryable"]


F = TypeVar("F", bound=Callable[..., Any])

DEFAULT_TIMEOUT = 10.0
DEFAULT_MAX_ATTEMPTS = 3


def _resolve_logger(args: tuple[Any, ...]) -> StructuredLogger:
    candidate = args[0] if args else None
    logger_attr = getattr(candidate, "logger", None)
    if isinstance(logger_attr, StructuredLogger):
        return logger_attr
    module = candidate.__class__.__module__ if candidate is not None else __name__
    return get_logger(module)


def retryable(
    *,
    max_attempts: int = DEFAULT_MAX_ATTEMPTS,
    timeout: float = DEFAULT_TIMEOUT,
    logger: StructuredLogger | None = None,
    stage: str = "http_request",
) -> Callable[[F], F]:
    """Retry ``func`` with exponential backoff on ``RequestException``."""

    def decorator(func: F) -> F:
        @functools.wraps(func)
        def call_with_timeout(*args: Any, **kwargs: Any) -> Any:
            if "timeout" not in kwargs:
                kwargs["timeout"] = timeout
            return func(*args, **kwargs)

        def _on_backoff(details: dict[str, Any]) -> None:
            log = logger or _resolve_logger(tuple(details.get("args", ())))
            exc = details["exception"]
            log.warning(
                "http_retry",
                stage=stage,
                attempt=details["tries"],
                wait=details["wait"],
                error=str(exc),
            )

        def _on_giveup(details: dict[str, Any]) -> None:
            log = logger or _resolve_logger(tuple(details.get("args", ())))
            exc = details["exception"]
            log.error(
                "http_fail",
                stage=stage,
                attempt=details["tries"],
                error=str(exc),
            )

        wrapped = backoff.on_exception(
            backoff.expo,
            requests.exceptions.RequestException,
            max_tries=max_attempts,
            on_backoff=_on_backoff,
            on_giveup=_on_giveup,
        )(call_with_timeout)

        return functools.wraps(func)(wrapped)  # type: ignore[return-value]

    return decorator
