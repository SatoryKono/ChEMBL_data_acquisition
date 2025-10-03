from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import IO, Iterator
import sys

from library.cli import LoggerConfig, configure_logger


class _TeeStream:
    """Mirror writes to multiple text streams."""

    def __init__(self, *streams: IO[str] | None) -> None:
        unique: list[IO[str]] = []
        for stream in streams:
            if stream is None:
                continue
            if any(stream is existing for existing in unique):
                continue
            unique.append(stream)
        self._streams: tuple[IO[str], ...] = tuple(unique)

    def write(self, data: str) -> int:  # pragma: no cover - trivial delegation
        for stream in self._streams:
            stream.write(data)
        return len(data)

    def flush(self) -> None:  # pragma: no cover - trivial delegation
        for stream in self._streams:
            flush = getattr(stream, "flush", None)
            if callable(flush):
                flush()

    def writable(self) -> bool:  # pragma: no cover - compatibility helper
        if not self._streams:
            return False
        return all(getattr(stream, "writable", lambda: True)() for stream in self._streams)

    def isatty(self) -> bool:  # pragma: no cover - compatibility helper
        return any(getattr(stream, "isatty", lambda: False)() for stream in self._streams)


@contextmanager
def configure_log_file(log_cfg: LoggerConfig, log_path: Path) -> Iterator[None]:
    """Configure ``log_cfg`` to tee structured logs to ``log_path`` and console."""

    log_path.parent.mkdir(parents=True, exist_ok=True)
    original_stream = log_cfg.stream
    with log_path.open("a", encoding="utf-8") as log_stream:
        tee_stream = _TeeStream(log_stream, original_stream)
        log_cfg.stream = tee_stream
        configure_logger(log_cfg)
        console_stream = original_stream or sys.stdout
        if console_stream is None:  # pragma: no cover - defensive fallback
            console_stream = sys.stdout
        print(
            f"[INFO] Structured logs are mirrored to '{log_path}'.",
            file=console_stream,
            flush=True,
        )
        try:
            yield
        finally:
            log_cfg.stream = original_stream
            configure_logger(log_cfg)


def build_log_path(script_name: str, *, directory: Path, timestamp: datetime | None = None) -> Path:
    """Return a predictable log path for ``script_name`` inside ``directory``."""

    moment = timestamp or datetime.now(timezone.utc)
    return directory / f"{script_name}_{moment.strftime('%Y%m%d')}.log"
