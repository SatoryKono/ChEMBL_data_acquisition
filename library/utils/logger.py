from __future__ import annotations

import logging
import sys
import time
from collections.abc import Iterable
from pathlib import Path
from typing import Any

__all__ = ["get_logger", "StructuredLogger"]


_LOG_FORMAT = "[%(asctime)s] [%(levelname)s] [%(name)s] %(message)s"
_LOG_DIR = Path("logs")


def _format_value(value: Any) -> str:
    if isinstance(value, str):
        escaped = value.replace("'", "\\'")
        return f"'{escaped}'"
    return repr(value)


def _build_message(
    event: str, message: str | None, kv: Iterable[tuple[str, Any]]
) -> str:
    parts = [message or event]
    for key, value in sorted(kv, key=lambda item: item[0]):
        parts.append(f"{key}={_format_value(value)}")
    return " ".join(parts)


class StructuredLogger:
    """Facade over :class:`logging.Logger` supporting structured calls."""

    def __init__(self, base_logger: logging.Logger) -> None:
        self._logger = base_logger

    @property
    def handlers(self) -> list[logging.Handler]:
        return list(self._logger.handlers)

    @property
    def name(self) -> str:
        return self._logger.name

    def addHandler(
        self, handler: logging.Handler
    ) -> None:  # pragma: no cover - delegation helper
        self._logger.addHandler(handler)

    def removeHandler(
        self, handler: logging.Handler
    ) -> None:  # pragma: no cover - delegation helper
        self._logger.removeHandler(handler)

    def _log(
        self,
        level: int,
        event: str,
        *args: Any,
        exc_info: Any | None = None,
        extra: dict[str, Any] | None = None,
        **kv: Any,
    ) -> None:
        if not self._logger.isEnabledFor(level):
            return

        message = event % args if args else event
        payload = dict(kv)
        if extra:
            payload.update(extra)
        formatted = _build_message(event, message, payload.items())
        self._logger.log(level, formatted, exc_info=exc_info)

    def log(self, level: int, event: str, *args: Any, **kwargs: Any) -> None:
        exc_info = kwargs.pop("exc_info", None)
        extra = kwargs.pop("extra", None)
        self._log(level, event, *args, exc_info=exc_info, extra=extra, **kwargs)

    def debug(self, event: str, *args: Any, **kwargs: Any) -> None:
        self.log(logging.DEBUG, event, *args, **kwargs)

    def info(self, event: str, *args: Any, **kwargs: Any) -> None:
        self.log(logging.INFO, event, *args, **kwargs)

    def warning(self, event: str, *args: Any, **kwargs: Any) -> None:
        self.log(logging.WARNING, event, *args, **kwargs)

    warn = warning

    def error(self, event: str, *args: Any, **kwargs: Any) -> None:
        self.log(logging.ERROR, event, *args, **kwargs)

    def exception(
        self, event: str, *args: Any, exc: BaseException | None = None, **kwargs: Any
    ) -> None:
        exc_info = kwargs.pop("exc_info", None)
        if exc_info is None:
            exc_info = exc
        self.log(logging.ERROR, event, *args, exc_info=exc_info or True, **kwargs)


def _ensure_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    return logger


def _configure_handlers(logger: logging.Logger, log_path: Path | None) -> None:
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        try:
            handler.close()
        except Exception:  # pragma: no cover - defensive close
            pass

    formatter = logging.Formatter(_LOG_FORMAT)
    formatter.converter = time.gmtime

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)


def get_logger(
    name: str, *, log_file: str | Path | None = None, level: int = logging.INFO
) -> StructuredLogger:
    """Return a structured logger writing to both stdout and ``log_file``."""

    logger = _ensure_logger(name)

    if log_file is None:
        log_path = _LOG_DIR / f"{name.replace('.', '_')}.log"
    else:
        log_path = Path(log_file)
        if log_path.is_dir():
            log_path = log_path / f"{name.replace('.', '_')}.log"

    _configure_handlers(logger, log_path)
    logger.setLevel(level)
    return StructuredLogger(logger)
