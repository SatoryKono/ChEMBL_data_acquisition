"""Shared command-line helpers."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any, Dict

from .config import Config, load_config


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


def build_parser(
    description: str, *, column: str, chunk_size: int = 10
) -> argparse.ArgumentParser:
    """Return an argument parser with shared options.

    Parameters
    ----------
    description:
        Text used in the parser description.
    column:
        Default column name for identifier extraction.
    chunk_size:
        Default chunk size for API requests.

    """
    parser = argparse.ArgumentParser(description=description)
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
    parser.add_argument(
        "--column",
        default=column,
        help="Identifier column in input CSV",
    )
    parser.add_argument("--sep", default=",", help="CSV delimiter")
    parser.add_argument("--encoding", default="utf8", help="File encoding")
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
    return parser


def configure_logging(
    level: str, *, fmt: str | None = None, datefmt: str | None = None
) -> None:
    """Configure root logging for command-line utilities.

    Parameters
    ----------
    level:
        Textual logging level (e.g. ``"INFO"``, ``"DEBUG"``).

    """
    numeric = getattr(logging, level.upper(), logging.INFO)
    logging.basicConfig(
        level=numeric,
        format=fmt or "%(levelname)s: %(message)s",
        datefmt=datefmt,
        force=True,
    )


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

    cfg = load_config(config_path, cli_overrides=cli_overrides)

    for arg, key in override_map.items():
        if not hasattr(args, arg):
            continue
        default = parser.get_default(arg)
        if getattr(args, arg) == default:
            setattr(args, arg, _get_cfg_value(cfg, key))

    return cfg


__all__ = ["build_parser", "configure_logging", "apply_config_overrides"]
