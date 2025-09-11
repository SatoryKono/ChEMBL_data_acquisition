"""Shared instance of :class:`~library.logging_setup.Logger`.

This module exposes a global logger configured with default settings. CLI
utilities should call :func:`library.cli.configure_logger` to update the
configuration (for example to set the log level or run identifier) before
emitting any log records.
"""

from __future__ import annotations

from .logging_setup import Logger, LoggerConfig, configure_logger

# Default logger used throughout the codebase.  The configuration is replaced
# by :func:`library.cli.configure_logger` when a CLI entry point starts up.
logger: Logger = configure_logger(LoggerConfig())

__all__ = ["logger"]
