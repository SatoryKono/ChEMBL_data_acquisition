"""Structured logging helpers for pipeline components."""

from __future__ import annotations

import logging
from contextlib import contextmanager
from contextvars import ContextVar
from time import perf_counter
from typing import Any, Iterator

_HANDLER_INSTALLED: bool = False

_EXTRA_CONTEXT: ContextVar[dict[str, Any]] = ContextVar("extra_context", default={})
_RUN_ID: ContextVar[str | None] = ContextVar("run_id", default=None)
_STAGE: ContextVar[str | None] = ContextVar("stage", default=None)
_STAGE_STARTED_AT: ContextVar[float | None] = ContextVar("stage_started_at", default=None)

_RESERVED_KWARGS = {"exc_info", "stack_info", "extra"}


def _now() -> float:
    """Return the current monotonic timestamp."""

    return perf_counter()


def _ensure_handler() -> None:
    """Install the structured formatter on the root logger once."""

    global _HANDLER_INSTALLED
    if _HANDLER_INSTALLED:
        return
    root = logging.getLogger()
    for handler in root.handlers:
        if isinstance(getattr(handler, "formatter", None), StructuredFormatter):
            _HANDLER_INSTALLED = True
            return
    handler = logging.StreamHandler()
    handler.setFormatter(StructuredFormatter())
    root.addHandler(handler)
    _HANDLER_INSTALLED = True


def _format_value(value: Any) -> str:
    """Render ``value`` for structured logging output."""

    if isinstance(value, float):
        return f"{value:.6f}".rstrip("0").rstrip(".") or "0"
    if isinstance(value, str):
        if value == "":
            return "''"
        if any(ch.isspace() for ch in value):
            escaped = value.replace("'", "\\'")
            return f"'{escaped}'"
        return value
    return repr(value)


class StructuredFormatter(logging.Formatter):
    """Formatter producing ``[LEVEL] [logger] message key=value`` records."""

    def format(self, record: logging.LogRecord) -> str:
        message = record.getMessage()
        head = f"[{record.levelname}] [{record.name}] {message}"
        data = getattr(record, "structured_data", {})
        if not isinstance(data, dict) or not data:
            return head
        parts = [head]
        for key in sorted(data):
            value = data[key]
            if value is None:
                continue
            parts.append(f"{key}={_format_value(value)}")
        return " ".join(parts)


class Logger(logging.LoggerAdapter):
    """Structured logger supporting contextual key-value data."""

    def __init__(self, logger: logging.Logger, extra: dict[str, Any] | None = None) -> None:
        super().__init__(logger, extra or {})

    def bind(self, **extra: Any) -> "Logger":
        """Return a new logger with ``extra`` context attached."""

        merged = {**self.extra, **extra}
        return Logger(self.logger, merged)

    def process(self, msg: Any, kwargs: dict[str, Any]) -> tuple[Any, dict[str, Any]]:
        structured: dict[str, Any] = {}
        if self.extra:
            structured.update(self.extra)
        provided_extra = kwargs.pop("extra", None)
        if isinstance(provided_extra, dict):
            structured.update(provided_extra.get("structured_data", {}))
            remaining = {k: v for k, v in provided_extra.items() if k != "structured_data"}
            structured.update(remaining)
        for key in list(kwargs.keys()):
            if key in _RESERVED_KWARGS:
                continue
            structured[key] = kwargs.pop(key)
        runtime_context = dict(_EXTRA_CONTEXT.get())
        if runtime_context:
            for key, value in runtime_context.items():
                structured.setdefault(key, value)
        run_id = _RUN_ID.get()
        if run_id is not None:
            structured.setdefault("run_id", run_id)
        stage = _STAGE.get()
        if stage is not None:
            structured.setdefault("stage", stage)
        stage_started = _STAGE_STARTED_AT.get()
        if stage_started is not None and "duration_s" not in structured:
            structured["duration_s"] = _now() - stage_started
        else:
            structured.setdefault("duration_s", None)
        kwargs["extra"] = {"structured_data": structured}
        return msg, kwargs

    @contextmanager
    def stage(self, name: str, **extra: Any) -> Iterator["Logger"]:
        """Bind ``stage`` for the lifetime of the context manager."""

        with log_context(stage=name, **extra):
            yield self


def get_logger(name: str) -> Logger:
    """Return a structured logger configured with the shared formatter."""

    _ensure_handler()
    return Logger(logging.getLogger(name))


@contextmanager
def log_context(
    *,
    run_id: str | None = None,
    stage: str | None = None,
    **extra: Any,
) -> Iterator[None]:
    """Temporarily bind contextual fields for structured logs."""

    tokens: list[tuple[ContextVar[Any], Any]] = []
    if extra:
        current = dict(_EXTRA_CONTEXT.get())
        current.update(extra)
        tokens.append((_EXTRA_CONTEXT, _EXTRA_CONTEXT.set(current)))
    if run_id is not None:
        tokens.append((_RUN_ID, _RUN_ID.set(run_id)))
    if stage is not None:
        tokens.append((_STAGE, _STAGE.set(stage)))
        tokens.append((_STAGE_STARTED_AT, _STAGE_STARTED_AT.set(_now())))
    try:
        yield
    finally:
        for var, token in reversed(tokens):
            var.reset(token)


@contextmanager
def log_stage(name: str, **extra: Any) -> Iterator[None]:
    """Shortcut for :func:`log_context` binding a stage name."""

    with log_context(stage=name, **extra):
        yield


__all__ = ["Logger", "StructuredFormatter", "get_logger", "log_context", "log_stage"]
