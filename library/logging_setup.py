"""Structured JSON logging utilities.

This module provides :class:`Logger` which emits structured log records as
JSON objects. Each log entry consists of a single line with a set of standard
fields and arbitrary key/value pairs supplied by the caller. Secrets contained
in the data are automatically redacted before output.
"""

from __future__ import annotations

import json
import logging
import sys
import threading
import time
import traceback
import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import IO, Any

# ``datetime.UTC`` was introduced in Python 3.11.
# Use ``timezone.utc`` for compatibility with earlier versions.
UTC = timezone.utc  # noqa: UP017

_SECRET_SUFFIXES: tuple[str, ...] = ("token", "key", "secret", "password")
_LEVELS = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARN": logging.WARNING,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
}

# Global lock ensuring thread-safe writes to the log stream.
_EMIT_LOCK = threading.Lock()

# Default attributes present on :class:`logging.LogRecord` instances.  Used to
# extract ``extra`` fields supplied via the standard :mod:`logging` API.
_LOG_RECORD_DEFAULT_KEYS = set(logging.LogRecord("", 0, "", 0, "", (), None).__dict__)


def _level_no(name: str) -> int:
    """Return numeric value for ``name``.

    Parameters
    ----------
    name:
        Level name such as ``"INFO"`` or ``"ERROR"``.  Unknown names default to
        ``logging.INFO``.
    """

    return _LEVELS.get(name.upper(), logging.INFO)


# ``slots`` is available from Python 3.10 onwards.  Supplying it on older
# versions raises ``TypeError``, so the argument is added conditionally.
_DATACLASS_KWARGS: dict[str, bool] = (
    {"slots": True} if sys.version_info >= (3, 10) else {}
)


@dataclass(**_DATACLASS_KWARGS)
class LoggerConfig:
    """Configuration for :class:`Logger`.

    Parameters
    ----------
    level:
        Minimum severity of emitted log records. Defaults to ``"INFO"``.
    run_id:
        Identifier for the current run or process. This value is included with
        every log record unless overridden via :meth:`Logger.bind`.
    redact_secrets:
        If ``True`` (default), values associated with keys ending in
        ``token``, ``key``, ``secret`` or ``password`` are replaced with
        ``"***"``.
    stream:
        Target text stream. Defaults to :data:`sys.stdout`.
    """

    level: str = "INFO"
    run_id: str = ""
    redact_secrets: bool = True
    stream: IO[str] = sys.stdout


class Logger:
    """Minimal structured logger writing JSON to a stream.

    Notes
    -----
    Instances are immutable. Methods that alter context such as
    :meth:`bind` return new logger objects.
    """

    def __init__(
        self, cfg: LoggerConfig, context: dict[str, Any] | None = None
    ) -> None:
        self._cfg = cfg
        self._context: dict[str, Any] = context or {}

    # ------------------------------------------------------------------
    # Internal helpers
    def _should_log(self, level: str) -> bool:
        return _level_no(level) >= _level_no(self._cfg.level)

    def _redact(self, data: dict[str, Any]) -> dict[str, Any]:
        if not self._cfg.redact_secrets:
            return data
        redacted = {}
        for k, v in data.items():
            if any(k.lower().endswith(s) for s in _SECRET_SUFFIXES):
                redacted[k] = "***"
            else:
                redacted[k] = v
        return redacted

    def _emit(self, record: dict[str, Any]) -> None:
        with _EMIT_LOCK:
            json.dump(record, self._cfg.stream)
            self._cfg.stream.write("\n")
            self._cfg.stream.flush()

    # ------------------------------------------------------------------
    # Public API
    def bind(self, **ctx: Any) -> Logger:
        """Return a new logger with ``ctx`` merged into the context."""

        new_ctx = {**self._context, **ctx}
        return Logger(self._cfg, new_ctx)

    def log(
        self,
        level: str,
        event: str,
        *,
        stage: str | None = None,
        msg: str | None = None,
        extra: dict[str, Any] | None = None,
        **kv: Any,
    ) -> None:
        """Emit a log record.

        Parameters
        ----------
        level:
            Severity of the event (e.g. ``"INFO"``).
        event:
            Short event name.
        stage:
            Optional stage name for pipelines. Overrides any bound ``stage``.
        msg:
            Optional human readable message.
        **kv:
            Additional key/value pairs added to the log record.
        """

        if not self._should_log(level):
            return

        record: dict[str, Any] = {
            "ts": datetime.now(UTC).isoformat(),
            "level": level.upper(),
            "event": event,
        }

        merged: dict[str, Any] = {**self._context, **kv}
        if extra:
            merged.update(extra)

        run_id = merged.pop("run_id", self._cfg.run_id)
        stage = stage if stage is not None else merged.pop("stage", None)

        record["run_id"] = run_id
        if stage is not None:
            record["stage"] = stage
        if msg is not None:
            record["msg"] = msg

        record.update(self._redact(merged))
        self._emit(record)

    # Convenience wrappers ------------------------------------------------
    def debug(
        self,
        event: str,
        *args: Any,
        extra: dict[str, Any] | None = None,
        **kv: Any,
    ) -> None:
        msg = event % args if args else None
        self.log("DEBUG", event, msg=msg, extra=extra, **kv)

    def info(
        self,
        event: str,
        *args: Any,
        extra: dict[str, Any] | None = None,
        **kv: Any,
    ) -> None:
        msg = event % args if args else None
        self.log("INFO", event, msg=msg, extra=extra, **kv)

    def warn(
        self,
        event: str,
        *args: Any,
        extra: dict[str, Any] | None = None,
        **kv: Any,
    ) -> None:
        msg = event % args if args else None
        self.log("WARN", event, msg=msg, extra=extra, **kv)

    # ``warning`` is an alias for compatibility with :mod:`logging` APIs.
    warning = warn

    def error(
        self,
        event: str,
        *args: Any,
        extra: dict[str, Any] | None = None,
        **kv: Any,
    ) -> None:
        msg = event % args if args else None
        self.log("ERROR", event, msg=msg, extra=extra, **kv)

    def exception(
        self,
        event: str,
        *args: Any,
        exc: BaseException | None = None,
        extra: dict[str, Any] | None = None,
        **kv: Any,
    ) -> None:
        """Log an exception at ``ERROR`` level.

        Parameters
        ----------
        event:
            Event name.
        exc:
            Exception instance. If ``None`` the current exception from
            :func:`sys.exc_info` is used.
        """

        if exc is None:
            exc = sys.exc_info()[1]
        if exc is None:  # pragma: no cover - defensive
            exc = Exception("unknown")
        msg = event % args if args else str(exc)
        tb = "".join(
            traceback.format_exception(exc.__class__, exc, exc.__traceback__, limit=3)
        ).strip()
        self.log(
            "ERROR",
            event,
            msg=msg,
            exc_type=exc.__class__.__name__,
            exc_message=str(exc),
            traceback=tb,
            extra=extra,
            **kv,
        )

    @contextmanager
    def stage(self, name: str, **kv: Any) -> Iterator[Logger]:
        """Context manager logging stage start/done/fail events.

        Parameters
        ----------
        name:
            Stage identifier used in event names and bound context.
        **kv:
            Extra context added for the duration of the stage.
        """

        bound = self.bind(stage=name, **kv)
        start = time.perf_counter()
        bound.info(f"{name}_start")
        try:
            yield bound
        except Exception as exc:  # pragma: no cover - branch exercised in tests
            bound.exception(f"{name}_fail", exc=exc)
            raise
        else:
            elapsed = time.perf_counter() - start
            bound.info(f"{name}_done", elapsed=elapsed)


# ---------------------------------------------------------------------------
# Public factory


_STRUCTURED_LOGGER_NAME = "library.structured"


def configure_logger(cfg: LoggerConfig, replace_root: bool = True) -> Logger:
    """Return a :class:`Logger` instance configured with ``cfg``.

    When ``replace_root`` is :data:`True`, the root :mod:`logging` logger is
    configured to forward records to the structured :class:`Logger` instance.
    This ensures a single log format is used for both direct :class:`Logger`
    calls and the standard :mod:`logging` module, including messages
    originating from :func:`warnings.warn` when :func:`logging.captureWarnings`
    is enabled. When ``replace_root`` is :data:`False`, the handler is instead
    attached to the ``"library.structured"`` logger leaving the root logger
    and ``py.warnings`` configuration untouched.

    Parameters
    ----------
    cfg:
        Logger configuration.
    replace_root:
        If :data:`True` (default), replace handlers on the root logger and the
        ``py.warnings`` logger. If :data:`False`, attach the handler to a
        dedicated ``"library.structured"`` logger.

    Returns
    -------
    Logger
        Configured logger instance.
    """

    logger = Logger(cfg)

    class _ForwardHandler(logging.Handler):
        """Forward ``logging`` records to :class:`Logger` as JSON lines."""

        def emit(self, record: logging.LogRecord) -> None:  # pragma: no cover - thin
            extras = {
                k: v
                for k, v in record.__dict__.items()
                if k not in _LOG_RECORD_DEFAULT_KEYS
            }
            extras["logger"] = record.name
            logger.log(record.levelname, record.getMessage(), extra=extras)

    handler = _ForwardHandler()
    if replace_root:
        root_logger = logging.getLogger()
        root_logger.handlers = [handler]
        root_logger.setLevel(_level_no(cfg.level))

        # Route ``warnings.warn`` calls through the structured logger.
        # ``logging`` redirects warnings to the ``py.warnings`` logger when
        # captureWarnings is enabled.  We attach the same handler used for the
        # root logger to ensure a single JSON-formatted output.
        logging.captureWarnings(True)
        warnings.simplefilter("default")
        warnings_logger = logging.getLogger("py.warnings")
        warnings_logger.handlers.clear()
        warnings_logger.addHandler(handler)
        warnings_logger.setLevel(_level_no(cfg.level))
        warnings_logger.propagate = False
    else:
        structured_logger = logging.getLogger(_STRUCTURED_LOGGER_NAME)
        structured_logger.handlers = [handler]
        structured_logger.setLevel(_level_no(cfg.level))
        structured_logger.propagate = False

    return logger
