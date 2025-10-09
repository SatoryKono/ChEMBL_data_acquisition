"""Helpers for configuring logging in CLI scripts."""

from __future__ import annotations

import logging
import os
import re
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import IO, Iterator

from ..common.logging_setup import LoggerConfig

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_DEFAULT_LOG_DIR = (_PROJECT_ROOT / "data" / "logs").resolve()

_SEPARATOR_PATTERN = re.compile(r"[\\/]+")
_SAFE_NAME_PATTERN = re.compile(r"[^0-9A-Za-z._-]+")


def _normalize_log_dir(path: Path | str) -> Path:
    """Return an absolute path for the requested log directory."""

    candidate = Path(path).expanduser()
    if candidate.is_absolute():
        return candidate.resolve()
    return (_PROJECT_ROOT / candidate).resolve()


def _default_log_dir() -> Path:
    """Return the default log directory considering environment overrides."""

    env_base = os.environ.get("CHEMBL_DA_BASE_PATH")
    if env_base:
        base_dir = _normalize_log_dir(env_base)
        return (base_dir / "logs").resolve()
    return _DEFAULT_LOG_DIR


def _current_date_str() -> str:
    """Return the current UTC date formatted as ``YYYYMMDD``."""

    return datetime.now(timezone.utc).strftime("%Y%m%d")


@dataclass(slots=True)
class CLILoggingContext:
    """Container returned by :func:`setup_cli_logging`."""

    log_path: Path
    log_cfg: LoggerConfig
    console_stream: IO[str]


def _normalise_script_name(script_name: str) -> str:
    """Return a stable identifier used as log file prefix."""

    candidate = script_name.strip()
    if not candidate:
        raise ValueError("script_name must not be empty")

    # Convert Windows-style separators before normalising through ``Path``.
    candidate = candidate.replace("\\", "/")
    name = Path(candidate).name

    if name.endswith(".py"):
        name = Path(name).stem

    # Replace unsupported characters to keep filenames portable across platforms.
    name = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    if not name or not any(ch.isalnum() for ch in name):
        return "pipeline"
    return name


@contextmanager
def setup_cli_logging(
    script_name: str | os.PathLike[str],
    log_cfg: LoggerConfig,
    date_str: str | None = None,
    *,
    log_dir: Path | None = None,
) -> Iterator[CLILoggingContext]:
    """Configure logging to mirror output to a file and the console."""

    if log_dir is not None:
        resolved_dir = _normalize_log_dir(log_dir)
    else:
        resolved_dir = _default_log_dir()
    resolved_dir.mkdir(parents=True, exist_ok=True)

    normalised_name = _normalise_script_name(script_name)

    if date_str:
        suffix = date_str
    else:
        suffix = _current_date_str()

    log_path = resolved_dir / f"{normalised_name}_{suffix}.log"
    log_path.touch(exist_ok=True)
    if not log_path.exists():  # pragma: no cover - defensive guard
        raise RuntimeError(f"Failed to create log file at '{log_path}'.")

    console_stream = getattr(log_cfg, "stream", None) or sys.stdout
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    handlers = list(log_cfg.handlers)
    handlers.append(file_handler)

    updated_cfg = LoggerConfig(
        level=log_cfg.level,
        run_id=log_cfg.run_id,
        generated_at=log_cfg.generated_at,
        redact_secrets=log_cfg.redact_secrets,
        stream=console_stream,
        handlers=handlers,
        logger_name=log_cfg.logger_name,
    )

    print(
        f"[INFO] Logs are written to '{log_path}'.",
        file=console_stream,
        flush=True,
    )

    try:
        yield CLILoggingContext(
            log_path=log_path,
            log_cfg=updated_cfg,
            console_stream=console_stream,
        )
    finally:
        file_handler.close()
