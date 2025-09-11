"""Shared command-line helpers.

Configuration loading errors are converted to user-facing messages using
``argparse.ArgumentParser.error``.
"""

from __future__ import annotations

import argparse
import uuid
from pathlib import Path
from typing import Any, Dict

from .config import Config, ConfigError, load_config
from . import log
from .logging_setup import Logger, LoggerConfig, configure_logger as _configure_logger


def create_logger_config(level: str) -> LoggerConfig:
    """Return :class:`LoggerConfig` with a random ``run_id``.

    Parameters
    ----------
    level:
        Desired logging level.

    Returns
    -------
    LoggerConfig
        Configuration containing ``run_id`` and ``level``.
    """

    return LoggerConfig(run_id=uuid.uuid4().hex, level=level)


def _positive_int(value: str) -> int:
    """Return ``value`` as a positive integer for ``argparse``.

    Parameters
    ----------
    value:
        String representation of the integer.

    Returns
    -------
    int
        The parsed positive integer.

    Raises
    ------
    argparse.ArgumentTypeError
        If ``value`` is not a positive integer.

    """
    try:
        ivalue = int(value)
    except ValueError as exc:  # pragma: no cover - handled by argparse
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("chunk size must be a positive integer")
    return ivalue


def add_common_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add shared CLI arguments to ``parser``.

    Parameters
    ----------
    parser:
        Parser to be extended with common arguments.

    Returns
    -------
    argparse.ArgumentParser
        The parser instance for convenience.
    """

    parser.add_argument("--log-level", default="INFO", help="Logging level")
    parser.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=Path("input.csv"),
        help="Input CSV file",
    )
    parser.add_argument(
        "--output",
        dest="output_csv",
        type=Path,
        default=None,
        help="Destination CSV file (default: auto-generate)",
    )
    parser.add_argument("--sep", default=",", help="CSV delimiter")
    parser.add_argument("--encoding", default="utf8", help="File encoding")
    return parser


def build_parser(
    description: str, *, column: str, chunk_size: int = 10
) -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Return a parser with shared options and logging configuration.

    Parameters
    ----------
    description:
        Text used in the parser description.
    column:
        Default column name for identifier extraction.
    chunk_size:
        Default chunk size for API requests.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        The parser and associated :class:`LoggerConfig`.
    """
    parser = argparse.ArgumentParser(description=description)
    add_common_arguments(parser)
    parser.add_argument(
        "--column",
        default=column,
        help="Identifier column in input CSV",
    )
    parser.add_argument(
        "--chunk-size",
        type=_positive_int,
        default=chunk_size,
        help="Maximum IDs per request",
    )
    parser.add_argument(
        "--config",
        dest="config",
        type=Path,
        default=Path("config.yaml"),
        help="YAML configuration file",
    )
    parser.add_argument(
        "--print-config",
        action="store_true",
        help="Print effective configuration and exit",
    )
    log_cfg = create_logger_config(parser.get_default("log_level"))
    return parser, log_cfg


def build_root_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Return a parser containing root-level options and logging config.

    The parser is created with ``add_help=False`` so it can be used as a parent
    for both the top-level parser and sub-commands, allowing shared options such
    as ``--config`` and ``--log-level`` to be supplied before or after the
    chosen sub-command.
    """

    parser = argparse.ArgumentParser(add_help=False)
    add_common_arguments(parser)
    parser.add_argument(
        "--config",
        dest="config",
        type=Path,
        default=Path("config.yaml"),
        help="YAML configuration file",
    )
    parser.add_argument(
        "--print-config",
        action="store_true",
        help="Print effective configuration and exit",
    )
    log_cfg = create_logger_config(parser.get_default("log_level"))
    return parser, log_cfg


def configure_logger(
    cfg: LoggerConfig, *, fmt: str | None = None, datefmt: str | None = None
) -> Logger:
    """Configure and return a structured logger based on ``cfg``.

    Parameters
    ----------
    cfg:
        Logging configuration containing ``run_id`` and ``level``.
    fmt, datefmt:
        Unused parameters retained for backward compatibility.

    Returns
    -------
    Logger
        Configured logger instance shared across the package.
    """

    # ``fmt`` and ``datefmt`` are ignored because JSON logs have a fixed
    # structure.  They are accepted to remain API compatible with the previous
    # implementation that configured :mod:`logging`.
    new_logger = _configure_logger(
        LoggerConfig(level=cfg.level, run_id=cfg.run_id, stream=cfg.stream)
    )
    # Update the shared logger instance in place so existing references remain
    # valid across the code base.
    log.logger._cfg = new_logger._cfg
    log.logger._context = {}
    return log.logger


# ---------------------------------------------------------------------------
# Configuration overrides
# ---------------------------------------------------------------------------

# Mapping of common CLI argument names to configuration paths.
_DEFAULT_OVERRIDES: Dict[str, str] = {
    "sep": "io.csv_sep",
    "encoding": "io.csv_encoding",
    "log_level": "log.level",
    "chunk_size": "jobs.chunk_size",
    "timeout": "api.timeout_read",
}


def _get_cfg_value(cfg: Config, path: str) -> Any:
    """Return the value in ``cfg`` located at ``path``.

    Parameters
    ----------
    cfg:
        Configuration object.
    path:
        Dot separated attribute path within ``cfg``.
    """

    current: Any = cfg
    for part in path.split("."):
        current = getattr(current, part)
    return current


def apply_config_overrides(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
    config_path: str | Path,
    mapping: Dict[str, str] | None = None,
) -> Config:
    """Load configuration applying command line overrides.

    This helper compares CLI arguments with parser defaults. Values that differ
    from the defaults are added to ``cli_overrides`` and passed to
    :func:`library.config.load_config`. After loading, ``args`` is updated with
    configuration values for options that were not explicitly provided.

    Parameters
    ----------
    args:
        Parsed command line arguments.
    parser:
        Argument parser used to determine default values.
    config_path:
        Location of the YAML configuration file.
    mapping:
        Optional mapping of argument names to ``Config`` attribute paths. The
        mapping is merged with a set of common defaults.

    Returns
    -------
    Config
        Loaded configuration object with overrides applied.

    Raises
    ------
    SystemExit
        If the configuration file cannot be loaded.
    """

    override_map = {**_DEFAULT_OVERRIDES, **(mapping or {})}

    cli_overrides: Dict[str, Any] = {}
    for arg, key in override_map.items():
        if not hasattr(args, arg):
            continue
        value = getattr(args, arg)
        default = parser.get_default(arg)
        if value != default:
            cli_overrides[key] = value

    try:
        cfg = load_config(config_path, cli_overrides=cli_overrides)
    except ConfigError as exc:
        log.logger.error("%s", exc)
        parser.error(str(exc))

    for arg, key in override_map.items():
        if not hasattr(args, arg):
            continue
        default = parser.get_default(arg)
        if getattr(args, arg) == default:
            setattr(args, arg, _get_cfg_value(cfg, key))

    return cfg


__all__ = [
    "LoggerConfig",
    "create_logger_config",
    "add_common_arguments",
    "build_parser",
    "build_root_parser",
    "configure_logger",
    "apply_config_overrides",
]
