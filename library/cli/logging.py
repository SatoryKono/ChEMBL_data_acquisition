"""Helpers for configuring structured logging in CLI scripts."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import IO, Iterator
import sys

from ..common.logging_setup import LoggerConfig

_DEFAULT_LOG_DIR = Path("data") / "logs"


class _TeeStream:
    """Mirror writes to multiple text streams without closing them."""

    def __init__(self, *streams: IO[str] | None) -> None:
        unique_streams: list[IO[str]] = []
        for stream in streams:
            if stream is None:
                continue
            if any(stream is existing for existing in unique_streams):
                continue
            unique_streams.append(stream)
        self._streams: tuple[IO[str], ...] = tuple(unique_streams)

    def write(self, data: str) -> int:  # pragma: no cover - direct delegation
        for stream in self._streams:
            stream.write(data)
        return len(data)

    def flush(self) -> None:  # pragma: no cover - direct delegation
        for stream in self._streams:
            flush = getattr(stream, "flush", None)
            if callable(flush):
                flush()

    def writable(self) -> bool:  # pragma: no cover - interface helper
        if not self._streams:
            return False
        return all(getattr(stream, "writable", lambda: True)() for stream in self._streams)

    def isatty(self) -> bool:  # pragma: no cover - interface helper
        return any(getattr(stream, "isatty", lambda: False)() for stream in self._streams)


def _current_date_str() -> str:
    """Return the current UTC date formatted as ``YYYYMMDD``."""

    return datetime.now(timezone.utc).strftime("%Y%m%d")


@dataclass(slots=True)
class CLILoggingContext:
    """Container returned by :func:`setup_cli_logging`."""

    log_path: Path
    log_cfg: LoggerConfig
    console_stream: IO[str]


@contextmanager
def setup_cli_logging(
    script_name: str,
    log_cfg: LoggerConfig,
    date_str: str | None = None,
    *,
    log_dir: Path | None = None,
) -> Iterator[CLILoggingContext]:
    """Configure logging to mirror structured output to a file and the console."""

    resolved_dir = Path(log_dir) if log_dir is not None else _DEFAULT_LOG_DIR
    resolved_dir.mkdir(parents=True, exist_ok=True)

    if date_str:
        suffix = date_str
    else:
        suffix = _current_date_str()

    log_path = resolved_dir / f"{script_name}_{suffix}.log"

    original_stream = getattr(log_cfg, "stream", None)
    with log_path.open("a", encoding="utf-8") as log_stream:
        tee_stream = _TeeStream(log_stream, original_stream)
        updated_cfg = LoggerConfig(
            level=log_cfg.level,
            run_id=log_cfg.run_id,
            redact_secrets=log_cfg.redact_secrets,
            stream=tee_stream,
        )
        console_stream = original_stream or sys.stdout
        if console_stream is None:  # pragma: no cover - defensive fallback
            console_stream = sys.stdout
        print(
            f"[INFO] Structured logs are mirrored to '{log_path}'.",
            file=console_stream,
            flush=True,
        )
        yield CLILoggingContext(
            log_path=log_path,
            log_cfg=updated_cfg,
            console_stream=console_stream,
        )
