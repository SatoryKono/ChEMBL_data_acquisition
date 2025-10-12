"""Command-line interface utilities."""

from .base import PipelineCLIBase
from .logging import setup_cli_logging
from .parser import (
    ConfigMetadata,
    Logger,
    LoggerConfig,
    add_common_arguments,
    apply_config_overrides,
    build_parser,
    build_root_parser,
    configure_logger,
    create_logger_config,
    path_argument,
    positive_int,
    prepare_io_paths,
    set_emit_legacy_help,
)

__all__ = [
    "Logger",
    "LoggerConfig",
    "add_common_arguments",
    "apply_config_overrides",
    "build_parser",
    "build_root_parser",
    "configure_logger",
    "create_logger_config",
    "set_emit_legacy_help",
    "PipelineCLIBase",
    "setup_cli_logging",
    "prepare_io_paths",
    "path_argument",
    "positive_int",
    "ConfigMetadata",
]
