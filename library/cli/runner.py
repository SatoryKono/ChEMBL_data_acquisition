"""Helpers for executing CLI commands with shared configuration scaffolding."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Mapping, Sequence
from importlib import import_module

from ..config import Config, ensure_dirs, print_config
from . import LoggerConfig, configure_logger


def run_cli_command(
    build_parser: Callable[[], tuple[argparse.ArgumentParser, LoggerConfig]],
    argv: Sequence[str] | None,
    config_mapping: Mapping[str, str] | None = None,
    handler: Callable[[Config, argparse.Namespace], int] | None = None,
) -> int:
    """Execute a CLI command using shared configuration and logging setup."""

    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)

    log_cfg.level = args.log_level
    pipeline_logger = configure_logger(log_cfg)
    pipeline_logger.info("pipeline_start", run_id=log_cfg.run_id)

    try:
        cli_module = import_module("library.cli")
        apply_overrides = getattr(cli_module, "apply_config_overrides")
        cfg = apply_overrides(
            args,
            parser,
            args.config,
            mapping=dict(config_mapping or {}),
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg)
            pipeline_logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        pipeline_logger = configure_logger(log_cfg)
    except (ValueError, TypeError) as exc:
        pipeline_logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        pipeline_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        pipeline_logger.error(
            "directory_setup_failed",
            error=str(exc),
        )
        pipeline_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    command = handler or getattr(args, "func", None)
    if command is None:
        pipeline_logger.error("missing_command_handler")
        pipeline_logger.info("pipeline_fail", run_id=log_cfg.run_id)
        raise ValueError("No CLI command handler available")

    exit_code = command(cfg, args)
    if exit_code == 0:
        pipeline_logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        pipeline_logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


__all__ = ["run_cli_command"]
