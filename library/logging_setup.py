"""Structured JSON logging utilities.

This module provides :class:`Logger` which emits structured log records as
JSON objects. Each log entry consists of a single line with a set of standard
fields and arbitrary key/value pairs supplied by the caller. Secrets contained
in the data are automatically redacted before output.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
import logging
import sys
import time
import traceback
from contextlib import contextmanager
from typing import Any, IO, Iterator, Optional


_SECRET_SUFFIXES: tuple[str, ...] = ("token", "key", "secret", "password")
_LEVELS = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARN": logging.WARNING,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
}


def _level_no(name: str) -> int:
    """Return numeric value for ``name``.

    Parameters
    ----------
    name:
        Level name such as ``"INFO"`` or ``"ERROR"``.  Unknown names default to
        ``logging.INFO``.
    """

    return _LEVELS.get(name.upper(), logging.INFO)


@dataclass(slots=True)
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
        self, cfg: LoggerConfig, context: Optional[dict[str, Any]] = None
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
        json.dump(record, self._cfg.stream)
        self._cfg.stream.write("\n")
        self._cfg.stream.flush()

    # ------------------------------------------------------------------
    # Public API
    def bind(self, **ctx: Any) -> "Logger":
        """Return a new logger with ``ctx`` merged into the context."""

        new_ctx = {**self._context, **ctx}
        return Logger(self._cfg, new_ctx)

    def log(
        self,
        level: str,
        event: str,
        *,
        stage: Optional[str] = None,
        msg: Optional[str] = None,
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
            "ts": datetime.now(timezone.utc).isoformat(),
            "level": level.upper(),
            "event": event,
        }

        merged: dict[str, Any] = {**self._context, **kv}

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
    def debug(self, event: str, **kv: Any) -> None:
        self.log("DEBUG", event, **kv)

    def info(self, event: str, **kv: Any) -> None:
        self.log("INFO", event, **kv)

    def warn(self, event: str, **kv: Any) -> None:
        self.log("WARN", event, **kv)

    def error(self, event: str, **kv: Any) -> None:
        self.log("ERROR", event, **kv)

    def exception(self, event: str, exc: BaseException, **kv: Any) -> None:
        """Log ``exc`` at ``ERROR`` level."""

        tb = "".join(
            traceback.format_exception(exc.__class__, exc, exc.__traceback__, limit=3)
        ).strip()
        self.log(
            "ERROR",
            event,
            msg=str(exc),
            exc_type=exc.__class__.__name__,
            exc_message=str(exc),
            traceback=tb,
            **kv,
        )

    @contextmanager
    def stage(self, name: str, **kv: Any) -> Iterator["Logger"]:
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
            bound.exception(f"{name}_fail", exc)
            raise
        else:
            elapsed = time.perf_counter() - start
            bound.info(f"{name}_done", elapsed=elapsed)


# ---------------------------------------------------------------------------
# Public factory


def configure_logger(cfg: LoggerConfig) -> Logger:
    """Return a :class:`Logger` instance configured with ``cfg``."""

    return Logger(cfg)
