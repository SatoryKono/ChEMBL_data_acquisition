"""Utilities for constructing command-line interfaces.

This module centralises shared CLI behaviour such as common argument
definitions, configuration loading and logging setup. Configuration errors are
presented to users via :meth:`argparse.ArgumentParser.error`.
"""

from __future__ import annotations

import argparse
import uuid
import os
import warnings
from datetime import datetime, timezone
from pathlib import Path, PureWindowsPath
from typing import Any

from pydantic import ValidationError

from .. import log
from ..config import Config, ConfigError, load_config
from ..logging_setup import Logger, LoggerConfig
from ..logging_setup import configure_logger as _configure_logger
from ..version import require_python_version
from ..utils.config import DEFAULT_CONFIG_PATH

require_python_version()


_DEPRECATED_OUTPUT_OPTION = "--output"
_DEPRECATED_OUTPUT_REPLACEMENT = "--final-out"
_DEPRECATED_OUTPUT_MESSAGE = (
    f"{_DEPRECATED_OUTPUT_OPTION} is deprecated; use {_DEPRECATED_OUTPUT_REPLACEMENT}"
)


def _emit_deprecated_output_warning() -> None:
    """Log a deprecation warning without failing when the logger stream is closed."""

    try:
        log.logger.warning(
            "cli_option_deprecated",
            option=_DEPRECATED_OUTPUT_OPTION,
            replacement=_DEPRECATED_OUTPUT_REPLACEMENT,
        )
    except ValueError as exc:
        if "closed file" not in str(exc).lower():
            raise
        warnings.warn(
            _DEPRECATED_OUTPUT_MESSAGE,
            DeprecationWarning,
            stacklevel=3,
        )


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


def positive_int(value: str) -> int:
    """Public wrapper around :func:`_positive_int` for CLI validators."""

    return _positive_int(value)


def _path_argument(value: str) -> Path:
    """Convert CLI ``value`` into a :class:`~pathlib.Path` instance.

    On non-Windows platforms command line inputs may still use Windows path
    separators (a backslash). This helper normalises such inputs by
    interpreting them as Windows paths and returning a POSIX-compatible
    ``Path``. The
    behaviour preserves regular POSIX paths and defers to :class:`Path` on
    Windows.
    """

    stripped = value.strip()
    if os.name != "nt" and "\\" in stripped:
        windows_path = PureWindowsPath(stripped)
        return Path(windows_path.as_posix())
    return Path(stripped)


def path_argument(value: str) -> Path:
    """Public wrapper around :func:`_path_argument` for CLI validators."""

    return _path_argument(value)


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
    When ``--final-out`` is omitted, a file named
    ``output.<input-stem>_<YYYYMMDD>.csv`` is created next to the input file.
    """

    log_level = "INFO" if defaults else argparse.SUPPRESS
    input_default: Path | object = Path("input.csv") if defaults else argparse.SUPPRESS
    output_default: Path | None | object = None if defaults else argparse.SUPPRESS
    final_default: Path | None | object = None if defaults else argparse.SUPPRESS
    sep_default: str | object = "," if defaults else argparse.SUPPRESS
    enc_default: str | object = "utf8" if defaults else argparse.SUPPRESS
    base_default: Path | None | object = None if defaults else argparse.SUPPRESS
    input_dir_default: Path | None | object = None if defaults else argparse.SUPPRESS
    output_dir_default: Path | None | object = None if defaults else argparse.SUPPRESS
    date_default: str | None | object = None if defaults else argparse.SUPPRESS
    force_default: bool | object = False if defaults else argparse.SUPPRESS
    skip_default: bool | object = False if defaults else argparse.SUPPRESS

    parser.add_argument("--log-level", default=log_level, help="Logging level")
    parser.add_argument(
        "--input",
        dest="input_csv",
        type=path_argument,
        default=input_default,
        help="Input CSV file",
    )
    parser.add_argument(
        "--final-out",
        dest="final_out",
        type=path_argument,
        default=final_default,
        help="Destination CSV file (default: output.<stem>_<YYYYMMDD>.csv)",
    )
    parser.add_argument(
        "--output",
        dest="output_csv",
        type=path_argument,
        default=output_default,
        help="Destination CSV file (default: output.<stem>_<YYYYMMDD>.csv)",
    )
    parser.add_argument(
        "--out",
        dest="_deprecated_out",
        type=path_argument,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--sep", default=sep_default, help="CSV delimiter")
    parser.add_argument("--encoding", default=enc_default, help="File encoding")
    parser.add_argument(
        "--base-path",
        dest="base_path",
        type=path_argument,
        default=base_default,
        help="Base directory for input and output data",
    )
    parser.add_argument(
        "--input-dir",
        dest="input_dir",
        type=path_argument,
        default=input_dir_default,
        help="Directory containing input artefacts",
    )
    parser.add_argument(
        "--output-dir",
        dest="output_dir",
        type=path_argument,
        default=output_dir_default,
        help="Directory receiving generated artefacts",
    )
    parser.add_argument(
        "--date",
        dest="date",
        default=date_default,
        help="Date prefix used when constructing default outputs",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=force_default,
        help="Overwrite outputs even when they already exist",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        default=skip_default,
        help="Skip processing if the destination file is present",
    )
    return parser


def build_parser(
    description: str,
    *,
    column: str,
    chunk_size: int = 10,
    size_option: str = "--chunk-size",
    size_dest: str = "chunk_size",
    size_help: str = "Maximum IDs per request",
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
    size_option:
        CLI option name controlling the batch size.
    size_dest:
        Attribute name assigned to the parsed batch size value.
    size_help:
        Help text shown for the batch size argument.

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
        size_option,
        dest=size_dest,
        type=_positive_int,
        default=chunk_size,
        help=size_help,
    )
    parser.add_argument(
        "--config",
        dest="config",
        type=path_argument,
        default=DEFAULT_CONFIG_PATH,
        help="YAML configuration file (default: config/config.yaml)",
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
        type=path_argument,
        default=DEFAULT_CONFIG_PATH,
        help="YAML configuration file (default: config/config.yaml)",
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
        type=path_argument,
        default=argparse.SUPPRESS,
        help="YAML configuration file (default: config/config.yaml)",
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
        Unsupported parameters retained for compatibility.

    Returns
    -------
    Logger
        Configured logger instance shared across the package.

    Raises
    ------
    ValueError
        If ``fmt`` or ``datefmt`` are supplied. Structured JSON logging has a
        fixed layout and cannot be customised.
    """

    if fmt is not None or datefmt is not None:
        raise ValueError(
            "Structured logging uses a fixed JSON layout; fmt/datefmt overrides"
            " are not supported."
        )

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
_PATH_PREFIXES: dict[str, list[str]] = {
    "api": ["sources", "chembl", "api"],
    "chembl": ["sources", "chembl", "cache"],
    "openalex": ["sources", "openalex"],
    "crossref": ["sources", "crossref"],
    "iuphar": ["sources", "iuphar"],
    "pubchem": ["sources", "pubchem"],
    "pubmed": ["sources", "pubmed"],
    "semantic_scholar": ["sources", "semantic_scholar"],
    "uniprot": ["sources", "uniprot", "api"],
    "uniprot_mapping": ["sources", "uniprot", "mapping"],
    "activity": ["sources", "chembl", "pipelines", "activity"],
    "assay": ["sources", "chembl", "pipelines", "assay"],
    "testitem": ["sources", "chembl", "pipelines", "testitem"],
    "document": ["sources", "chembl", "pipelines", "document"],
    "target": ["sources", "chembl", "pipelines", "target"],
    "io": ["local", "io"],
    "resources": ["local", "resources"],
    "init": ["local", "init"],
    "log": ["system", "log"],
    "rate": ["system", "rate"],
    "retry": ["system", "retry"],
    "doc_type": ["system", "doc_type"],
}


def _normalize_path(path: str) -> str:
    parts = path.split(".") if path else []
    if not parts:
        return path
    head, tail = parts[0], parts[1:]
    prefix = _PATH_PREFIXES.get(head)
    if not prefix:
        return path
    new_parts = prefix + tail
    return ".".join(new_parts)


_DEFAULT_OVERRIDES: dict[str, str] = {
    key: _normalize_path(value)
    for key, value in {
        "base_path": "io.base_path",
        "input_dir": "io.input_dir",
        "output_dir": "io.output_dir",
        "cache_dir": "io.cache_dir",
        "sep": "io.csv_sep",
        "encoding": "io.csv_encoding",
        "log_level": "log.level",
        "timeout": "api.timeout_read",
    }.items()
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

    override_map = {
        arg: _normalize_path(path)
        for arg, path in {**_DEFAULT_OVERRIDES, **(mapping or {})}.items()
    }

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
        cfg = load_config(
            config_path,
            cli_overrides=cli_overrides,
            strict=True,
        )
    except ConfigError as exc:
        log.logger.error(
            "config_load_failed",
            error=str(exc),
            config=str(config_path),
        )
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


def _resolve_base_path(value: Path | str | None) -> Path | None:
    if value is None or value is argparse.SUPPRESS:
        return None
    path = (value if isinstance(value, Path) else Path(value)).expanduser()
    if path.is_absolute():
        return path
    return (Path.cwd() / path).resolve()


def _resolve_directory(
    directory: Path | str | None,
    *,
    base: Path | None,
) -> Path | None:
    if directory is None or directory is argparse.SUPPRESS:
        return None
    resolved = (
        directory if isinstance(directory, Path) else Path(directory)
    ).expanduser()
    if resolved.is_absolute():
        return resolved
    if base is not None:
        return (base / resolved).resolve()
    return resolved.resolve()


def _resolve_file(
    path: Path | str | None,
    *,
    directory: Path | None,
    base: Path | None,
) -> Path | None:
    if path is None or path is argparse.SUPPRESS:
        return None
    candidate = (path if isinstance(path, Path) else Path(path)).expanduser()
    if candidate.is_absolute():
        return candidate
    if directory is not None:
        return (directory / candidate).resolve()
    if base is not None:
        return (base / candidate).resolve()
    return candidate


def _raw_suffix(raw_format: str, fallback: str) -> str:
    value = raw_format.strip().lower()
    if not value:
        return fallback
    if value == "csv":
        return fallback
    if value == "parquet":
        return ".parquet"
    if value.startswith("."):
        return value
    return f".{value}"


def prepare_io_paths(
    args: argparse.Namespace,
    *,
    input_default: str | Path | None = None,
    output_stem: str | None = None,
    suffix: str = ".csv",
) -> None:
    """Normalize CLI namespace paths and populate derived locations."""

    base_path = _resolve_base_path(getattr(args, "base_path", None))
    setattr(args, "base_path", base_path)

    input_dir = _resolve_directory(getattr(args, "input_dir", None), base=base_path)
    setattr(args, "input_dir", input_dir)

    output_dir = _resolve_directory(getattr(args, "output_dir", None), base=base_path)
    setattr(args, "output_dir", output_dir)

    cache_dir = _resolve_directory(getattr(args, "cache_dir", None), base=base_path)
    setattr(args, "cache_dir", cache_dir)

    deprecated_out = getattr(args, "_deprecated_out", argparse.SUPPRESS)
    if deprecated_out is not argparse.SUPPRESS:
        current_output = getattr(args, "output_csv", None)
        if current_output in (None, argparse.SUPPRESS):
            setattr(args, "output_csv", deprecated_out)
        setattr(args, "_out_alias_used", True)

    current_input = getattr(args, "input_csv", None)
    if current_input in (None, argparse.SUPPRESS) and input_default is not None:
        current_input = Path(input_default)
    resolved_input = _resolve_file(current_input, directory=input_dir, base=base_path)
    if resolved_input is not None:
        setattr(args, "input_csv", resolved_input)

    final_candidate = getattr(args, "final_out", None)
    output_candidate = getattr(args, "output_csv", None)

    if output_candidate not in (None, argparse.SUPPRESS):
        if final_candidate in (None, argparse.SUPPRESS):
            final_candidate = output_candidate
            setattr(args, "final_out", final_candidate)
        _emit_deprecated_output_warning()

    raw_format_value = getattr(args, "raw_format", None)
    if raw_format_value in (None, argparse.SUPPRESS):
        raw_format_str = "csv"
    else:
        raw_format_str = str(raw_format_value).lower()
    setattr(args, "raw_format", raw_format_str)

    final_value = getattr(args, "final_out", None)
    if final_value in (None, argparse.SUPPRESS):
        final_value = None

    resolved_output = _resolve_file(
        final_value,
        directory=output_dir,
        base=base_path,
    )

    date_value = getattr(args, "date", None)
    if isinstance(date_value, str):
        date_str: str | None = date_value
    else:
        date_str = None

    if resolved_output is None and output_stem is not None:
        target_dir = output_dir or base_path
        if target_dir is None and resolved_input is not None:
            target_dir = resolved_input.parent
        if target_dir is not None:
            effective_date = date_str or datetime.now(timezone.utc).strftime("%Y%m%d")
            filename = f"output.{output_stem}_{effective_date}{suffix}"
            resolved_output = (target_dir / filename).resolve()
            date_str = effective_date

    if resolved_output is not None:
        resolved_output = _resolve_file(
            resolved_output,
            directory=output_dir,
            base=base_path,
        )
        setattr(args, "output_csv", resolved_output)
        setattr(args, "final_out", resolved_output)
    elif isinstance(final_value, (str, Path)):
        candidate_path = Path(final_value)
        setattr(args, "output_csv", candidate_path)
        setattr(args, "final_out", candidate_path)

    raw_value = getattr(args, "raw_out", None)
    if raw_value in (None, argparse.SUPPRESS):
        if resolved_output is not None:
            if raw_format_str == "csv":
                resolved_raw = resolved_output
            else:
                suffix_value = _raw_suffix(raw_format_str, suffix)
                resolved_raw = resolved_output.with_suffix(suffix_value)
        else:
            resolved_raw = None
    else:
        resolved_raw = _resolve_file(raw_value, directory=output_dir, base=base_path)

    if resolved_raw is not None:
        resolved_raw = _resolve_file(
            resolved_raw,
            directory=output_dir,
            base=base_path,
        )
    setattr(args, "raw_out", resolved_raw)

    if date_str is not None:
        setattr(args, "date", date_str)

    setattr(args, "force", bool(getattr(args, "force", False)))
    setattr(args, "skip_existing", bool(getattr(args, "skip_existing", False)))


__all__ = [
    "LoggerConfig",
    "Logger",
    "create_logger_config",
    "add_common_arguments",
    "build_parser",
    "build_root_parser",
    "configure_logger",
    "apply_config_overrides",
    "prepare_io_paths",
]
