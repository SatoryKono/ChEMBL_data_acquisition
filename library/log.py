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
#
# The logger is bound with ``status`` and ``rps`` fields to ensure that each
# emitted record contains these keys. Callers can override them by supplying
# values via the ``extra`` argument when logging.
logger: Logger = configure_logger(LoggerConfig()).bind(status=None, rps=None)

__all__ = ["logger"]
