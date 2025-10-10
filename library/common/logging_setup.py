"""Utilities for configuring standard logging across CLI tools."""

from __future__ import annotations

import logging
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import IO, Any

from .run_context import RunContext, set_current

_LEVELS: dict[str, int] = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARN": logging.WARNING,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
}

_SECRET_SUFFIXES: tuple[str, ...] = ("token", "key", "secret", "password")

_FORMAT = "[%(asctime)s] [%(levelname)s] [%(name)s] %(message)s"


def _close_handlers(handlers: list[logging.Handler]) -> None:
    for handler in handlers:
        try:
            handler.close()
        except Exception:  # pragma: no cover - defensive
            pass


def _level_no(name: str) -> int:
    """Return numeric logging level for ``name``.

    Unknown names default to :data:`logging.INFO`.
    """

    return _LEVELS.get(name.upper(), logging.INFO)


@dataclass(slots=True)
class LoggerConfig:
    """Configuration shared by :class:`Logger` and :func:`configure_logger`."""

    level: str = "INFO"
    run_id: str = ""
    generated_at: str = ""
    redact_secrets: bool = True
    stream: IO[str] | None = sys.stdout
    handlers: list[logging.Handler] = field(default_factory=list)
    logger_name: str = "chembl"


def _format_value(value: Any) -> str:
    if isinstance(value, str):
        escaped = value.replace("'", "\\'")
        return f"'{escaped}'"
    return repr(value)


def _format_message(event: str, message: str | None, context: dict[str, Any]) -> str:
    parts: list[str] = [event]
    if message and message != event:
        parts.append(f"- {message}")
    for key in sorted(context):
        parts.append(f"{key}={_format_value(context[key])}")
    return " ".join(parts)


class Logger:
    """Thin wrapper over :mod:`logging` providing structured helpers."""

    def __init__(
        self,
        cfg: LoggerConfig,
        *,
        base_logger: logging.Logger | None = None,
        context: dict[str, Any] | None = None,
    ) -> None:
        self._cfg = cfg
        self._logger = base_logger or logging.getLogger(cfg.logger_name)
        self._context: dict[str, Any] = context.copy() if context else {}

    def bind(self, **ctx: Any) -> Logger:
        merged = {**self._context, **ctx}
        return Logger(self._cfg, base_logger=self._logger, context=merged)

    def _redact(self, data: dict[str, Any]) -> dict[str, Any]:
        if not self._cfg.redact_secrets:
            return data
        redacted: dict[str, Any] = {}
        for key, value in data.items():
            if any(key.lower().endswith(suffix) for suffix in _SECRET_SUFFIXES):
                redacted[key] = "***"
            else:
                redacted[key] = value
        return redacted

    def log(
        self,
        level: str,
        event: str,
        *,
        stage: str | None = None,
        msg: str | None = None,
        extra: dict[str, Any] | None = None,
        exc_info: Any | None = None,
        **kv: Any,
    ) -> None:
        level_no = _level_no(level)
        if not self._logger.isEnabledFor(level_no):
            return

        context = dict(self._context)
        if stage is not None:
            context["stage"] = stage
        if extra:
            context.update(extra)
        context.update(kv)
        if self._cfg.run_id and "run_id" not in context:
            context["run_id"] = self._cfg.run_id

        context = self._redact(context)
        message = _format_message(event, msg, context)

        log_kwargs: dict[str, Any] = {}
        if exc_info is not None:
            if exc_info is True:
                log_kwargs["exc_info"] = True
            elif isinstance(exc_info, BaseException):
                log_kwargs["exc_info"] = (
                    exc_info.__class__,
                    exc_info,
                    exc_info.__traceback__,
                )
            else:
                log_kwargs["exc_info"] = exc_info

        self._logger.log(level_no, message, **log_kwargs)

    def debug(self, event: str, *args: Any, **kv: Any) -> None:
        msg = event % args if args else None
        exc_info = kv.pop("exc_info", None)
        self.log("DEBUG", event, msg=msg, exc_info=exc_info, **kv)

    def info(self, event: str, *args: Any, **kv: Any) -> None:
        msg = event % args if args else None
        exc_info = kv.pop("exc_info", None)
        self.log("INFO", event, msg=msg, exc_info=exc_info, **kv)

    def warn(self, event: str, *args: Any, **kv: Any) -> None:
        msg = event % args if args else None
        exc_info = kv.pop("exc_info", None)
        self.log("WARN", event, msg=msg, exc_info=exc_info, **kv)

    warning = warn

    def error(self, event: str, *args: Any, **kv: Any) -> None:
        msg = event % args if args else None
        exc_info = kv.pop("exc_info", None)
        self.log("ERROR", event, msg=msg, exc_info=exc_info, **kv)

    def exception(
        self,
        event: str,
        *args: Any,
        exc: BaseException | None = None,
        extra: dict[str, Any] | None = None,
        **kv: Any,
    ) -> None:
        if exc is None:
            exc = sys.exc_info()[1]
        if exc is None:
            exc = Exception(event)
        msg = event % args if args else str(exc)
        context = {"exception": exc.__class__.__name__, "error": str(exc)}
        if extra:
            context.update(extra)
        else:
            extra = context
        self.log(
            "ERROR",
            event,
            msg=msg,
            extra=context,
            exc_info=(exc.__class__, exc, exc.__traceback__),
            **kv,
        )

    @contextmanager
    def stage(self, name: str, **kv: Any) -> Iterator[Logger]:
        bound = self.bind(stage=name, **kv)
        start = time.perf_counter()
        bound.info(f"{name}_start")
        try:
            yield bound
        except Exception as exc:  # pragma: no cover - defensive
            bound.exception(f"{name}_fail", exc=exc)
            raise
        else:
            elapsed = time.perf_counter() - start
            bound.info(f"{name}_done", elapsed=elapsed)


def configure_logger(cfg: LoggerConfig, replace_root: bool = True) -> Logger:
    """Configure the logging subsystem and return a :class:`Logger`."""

    level_no = _level_no(cfg.level)
    formatter = logging.Formatter(_FORMAT)
    # Ensure timestamps are rendered in UTC to keep logs deterministic across environments
    formatter.converter = time.gmtime

    handlers: list[logging.Handler] = []
    if cfg.stream is not None:
        handlers.append(logging.StreamHandler(cfg.stream))
    handlers.extend(cfg.handlers)

    for handler in handlers:
        handler.setLevel(level_no)
        handler.setFormatter(formatter)

    if replace_root:
        root_logger = logging.getLogger()
        _close_handlers(list(root_logger.handlers))
        root_logger.handlers = []
        for handler in handlers:
            root_logger.addHandler(handler)
        root_logger.setLevel(level_no)
        logging.captureWarnings(True)
        warnings_logger = logging.getLogger("py.warnings")
        _close_handlers(list(warnings_logger.handlers))
        warnings_logger.handlers = []
        for handler in handlers:
            warnings_logger.addHandler(handler)
        warnings_logger.setLevel(level_no)
        warnings_logger.propagate = False
    else:
        target = logging.getLogger(cfg.logger_name)
        _close_handlers(list(target.handlers))
        target.handlers = []
        for handler in handlers:
            target.addHandler(handler)
        target.setLevel(level_no)

    base_logger = logging.getLogger(cfg.logger_name)
    base_logger.setLevel(level_no)
    base_logger.propagate = True

    set_current(RunContext(run_id=cfg.run_id, generated_at=cfg.generated_at))
    return Logger(cfg, base_logger=base_logger)


__all__ = ["Logger", "LoggerConfig", "configure_logger"]
