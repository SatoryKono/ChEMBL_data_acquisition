"""Shared command-line helpers.

The parsers built by :func:`build_parser` expose a ``--config`` option whose
value can be passed to :func:`library.config.load_config`.  If the configuration
file is missing or malformed, calling :func:`~library.config.load_config` will
raise :class:`library.config.ConfigError`.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .config import DEFAULT_PATH


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

    Notes
    -----
    The returned parser defines common arguments used across scripts,
    including ``--config`` for specifying the path to a YAML configuration
    file.  Consumers should load the file with
    :func:`library.config.load_config` and handle
    :class:`library.config.ConfigError` for missing or malformed files.

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
        type=Path,
        default=DEFAULT_PATH,
        help="Path to configuration YAML file",
    )
    return parser


def configure_logging(level: str) -> None:
    """Configure root logging for command-line utilities.

    Parameters
    ----------
    level:
        Textual logging level (e.g. ``"INFO"``, ``"DEBUG"``).

    """
    numeric = getattr(logging, level.upper(), logging.INFO)
    logging.basicConfig(level=numeric, force=True)
