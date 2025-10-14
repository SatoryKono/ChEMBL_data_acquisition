"""Structured logging utilities with contextual metadata support."""

from __future__ import annotations

import logging
import sys
import time
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any

from contextvars import ContextVar

__all__ = ["StructuredLogger", "get_logger", "log_context"]


_LOG_FORMAT = "[%(levelname)s] [%(name)s] %(message)s"


def _format_value(value: Any) -> str:
    if isinstance(value, str):
        escaped = value.replace("'", "\\'")
        return f"'{escaped}'"
    return repr(value)


def _build_message(message: str, kv: Iterable[tuple[str, Any]]) -> str:
    parts = [message]
    for key, value in sorted(kv, key=lambda item: item[0]):
        parts.append(f"{key}={_format_value(value)}")
    return " ".join(parts)


@dataclass(slots=True)
class _LogContext:
    run_id: str | None = None
    stage: str | None = None
    started_at: float | None = None


_CONTEXT: ContextVar[_LogContext] = ContextVar(
    "structured_logging_context", default=_LogContext()
)


def _resolve_context() -> _LogContext:
    ctx = _CONTEXT.get()
    # ``ContextVar`` returns the same default instance for every access, so
    # clone the context to avoid accidental mutation in nested scopes.
    return _LogContext(run_id=ctx.run_id, stage=ctx.stage, started_at=ctx.started_at)


class StructuredLogger:
    """Thin facade over :mod:`logging` supporting structured output."""

    def __init__(self, base_logger: logging.Logger) -> None:
        self._logger = base_logger

    @property
    def handlers(self) -> list[logging.Handler]:
        return list(self._logger.handlers)

    def addHandler(self, handler: logging.Handler) -> None:  # pragma: no cover - delegation
        self._logger.addHandler(handler)

    def removeHandler(
        self, handler: logging.Handler
    ) -> None:  # pragma: no cover - delegation
        self._logger.removeHandler(handler)

    def _log(
        self,
        level: int,
        message: str,
        *args: Any,
        exc_info: Any | None = None,
        extra: dict[str, Any] | None = None,
        **kv: Any,
    ) -> None:
        if not self._logger.isEnabledFor(level):
            return

        payload = dict(kv)
        if extra:
            payload.update(extra)

        context = _resolve_context()
        if context.run_id is not None and "run_id" not in payload:
            payload["run_id"] = context.run_id
        if context.stage is not None and "stage" not in payload:
            payload["stage"] = context.stage
        if context.started_at is not None and "duration_s" not in payload:
            payload["duration_s"] = round(time.perf_counter() - context.started_at, 6)

        formatted = _build_message(message % args if args else message, payload.items())
        self._logger.log(level, formatted, exc_info=exc_info)

    def debug(self, message: str, *args: Any, **kwargs: Any) -> None:
        self._log(logging.DEBUG, message, *args, **kwargs)

    def info(self, message: str, *args: Any, **kwargs: Any) -> None:
        self._log(logging.INFO, message, *args, **kwargs)

    def warning(self, message: str, *args: Any, **kwargs: Any) -> None:
        self._log(logging.WARNING, message, *args, **kwargs)

    warn = warning

    def error(self, message: str, *args: Any, **kwargs: Any) -> None:
        self._log(logging.ERROR, message, *args, **kwargs)

    def exception(
        self,
        message: str,
        *args: Any,
        exc: BaseException | None = None,
        **kwargs: Any,
    ) -> None:
        exc_info = kwargs.pop("exc_info", None)
        if exc_info is None:
            exc_info = exc or True
        self._log(logging.ERROR, message, *args, exc_info=exc_info, **kwargs)

    @contextmanager
    def context(
        self, *, run_id: str | None = None, stage: str | None = None
    ) -> Iterator[StructuredLogger]:
        with log_context(run_id=run_id, stage=stage):
            yield self


def get_logger(name: str) -> StructuredLogger:
    """Return a structured logger configured for stdout output."""

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(_LOG_FORMAT))
        logger.addHandler(handler)

    return StructuredLogger(logger)


@contextmanager
def log_context(
    *, run_id: str | None = None, stage: str | None = None
) -> Iterator[None]:
    """Temporarily bind contextual metadata for structured logs."""

    current = _resolve_context()
    next_context = _LogContext(
        run_id=run_id if run_id is not None else current.run_id,
        stage=stage if stage is not None else current.stage,
        started_at=time.perf_counter(),
    )
    token = _CONTEXT.set(next_context)
    try:
        yield
    finally:
        _CONTEXT.reset(token)
