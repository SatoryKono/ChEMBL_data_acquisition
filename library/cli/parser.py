"""Utilities for constructing command-line interfaces.

This module centralises shared CLI behaviour such as common argument
definitions, configuration loading and logging setup. Configuration errors are
presented to users via :meth:`argparse.ArgumentParser.error`.
"""

from __future__ import annotations

import argparse
import uuid
from pathlib import Path
from typing import Any

from pydantic import ValidationError

from .. import log
from ..config import Config, ConfigError, load_config
from ..logging_setup import Logger, LoggerConfig
from ..logging_setup import configure_logger as _configure_logger
from ..version import require_python_version

require_python_version()


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


def add_common_arguments(
    parser: argparse.ArgumentParser, *, defaults: bool = True
) -> argparse.ArgumentParser:
    """Add shared CLI arguments to ``parser``.

    Parameters
    ----------
    parser:
        Parser to be extended with common arguments.
    defaults:
        Whether to apply default values. When ``False`` the options are added
        with ``argparse.SUPPRESS`` so that they do not override values provided
        on a parent parser. This is useful for sub-commands that share global
        options while still allowing those options to be specified before the
        sub-command name.

    Returns
    -------
    argparse.ArgumentParser
        The parser instance for convenience.

    Notes
    -----
    When ``--output`` is omitted, a file named
    ``output_<input-stem>_<YYYYMMDD>.csv`` is created next to the input file.
    """

    log_level = "INFO" if defaults else argparse.SUPPRESS
    input_default: Path | object = Path("input.csv") if defaults else argparse.SUPPRESS
    output_default: Path | None | object = None if defaults else argparse.SUPPRESS
    sep_default: str | object = "," if defaults else argparse.SUPPRESS
    enc_default: str | object = "utf8" if defaults else argparse.SUPPRESS

    parser.add_argument("--log-level", default=log_level, help="Logging level")
    parser.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=input_default,
        help="Input CSV file",
    )
    parser.add_argument(
        "--output",
        dest="output_csv",
        type=Path,
        default=output_default,
        help="Destination CSV file (default: output_<stem>_<YYYYMMDD>.csv)",
    )
    parser.add_argument("--sep", default=sep_default, help="CSV delimiter")
    parser.add_argument("--encoding", default=enc_default, help="File encoding")
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


def build_root_parser() -> (
    tuple[argparse.ArgumentParser, argparse.ArgumentParser, LoggerConfig]
):
    """Return parsers containing root-level options and logging config.

    Two parsers are produced:

    ``root``
        Includes default values and is intended as the parent for the top-level
        parser.
    ``shared``
        Contains the same arguments but without defaults so that sub-commands
        can inherit the options without overriding values supplied before the
        sub-command name.

    Returns
    -------
    tuple[argparse.ArgumentParser, argparse.ArgumentParser, LoggerConfig]
        The parser with defaults, the version without defaults for sub-commands
        and the associated :class:`LoggerConfig`.
    """

    root = argparse.ArgumentParser(add_help=False)
    add_common_arguments(root)
    root.add_argument(
        "--config",
        dest="config",
        type=Path,
        default=Path("config.yaml"),
        help="YAML configuration file",
    )
    root.add_argument(
        "--print-config",
        action="store_true",
        help="Print effective configuration and exit",
    )

    shared = argparse.ArgumentParser(add_help=False)
    add_common_arguments(shared, defaults=False)
    shared.add_argument(
        "--config",
        dest="config",
        type=Path,
        default=argparse.SUPPRESS,
        help="YAML configuration file",
    )
    shared.add_argument(
        "--print-config",
        action="store_true",
        help="Print effective configuration and exit",
    )

    log_cfg = create_logger_config(root.get_default("log_level"))
    return root, shared, log_cfg


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
    # structure. They are accepted to remain API compatible with the previous
    # implementation that configured :mod:`logging`.
    new_logger = _configure_logger(
        LoggerConfig(level=cfg.level, run_id=cfg.run_id, stream=cfg.stream)
    )
    # Update the shared logger instance in place so existing references remain
    # valid across the code base.
    log.logger._cfg = new_logger._cfg
    # Preserve default structured fields for downstream log records.
    log.logger._context = {"status": None, "rps": None}
    return log.logger


# ---------------------------------------------------------------------------
# Configuration overrides
# ---------------------------------------------------------------------------

# Mapping of common CLI argument names to configuration paths.
_DEFAULT_OVERRIDES: dict[str, str] = {
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
    mapping: dict[str, str] | None = None,
    *,
    base_parser: argparse.ArgumentParser | None = None,
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
        Argument parser used to determine default values for command specific
        options.
    config_path:
        Location of the YAML configuration file.
    mapping:
        Optional mapping of argument names to ``Config`` attribute paths. The
        mapping is merged with a set of common defaults.
    base_parser:
        Optional parser providing fallback defaults for shared arguments. This
        is useful when ``parser`` is a sub-command parser without default values
        for global options.

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

    cli_overrides: dict[str, Any] = {}
    for arg, key in override_map.items():
        if not hasattr(args, arg):
            continue
        value = getattr(args, arg)
        default = parser.get_default(arg)
        if base_parser is not None and (
            default is None or default is argparse.SUPPRESS
        ):
            default = base_parser.get_default(arg)
        if value != default:
            cli_overrides[key] = value

    try:
        cfg = load_config(config_path, cli_overrides=cli_overrides)
    except ConfigError as exc:
        log.logger.error("%s", exc)
        parser.error(str(exc))
    except ValidationError as exc:
        raise ValueError(str(exc)) from exc

    for arg, key in override_map.items():
        if not hasattr(args, arg):
            continue
        default = parser.get_default(arg)
        if base_parser is not None and (
            default is None or default is argparse.SUPPRESS
        ):
            default = base_parser.get_default(arg)
        if getattr(args, arg) == default:
            setattr(args, arg, _get_cfg_value(cfg, key))

    return cfg


__all__ = [
    "LoggerConfig",
    "Logger",
    "create_logger_config",
    "add_common_arguments",
    "build_parser",
    "build_root_parser",
    "configure_logger",
    "apply_config_overrides",
]
