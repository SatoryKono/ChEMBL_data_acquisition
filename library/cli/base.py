"""Base classes for command-line pipeline entry points."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .logging import CLILoggingContext, setup_cli_logging
from .parser import Logger, LoggerConfig, configure_logger
from ..common.log import logger as default_logger


class PipelineCLIBase:
    """Template for command-line interfaces orchestrating data pipelines."""

    def build_parser(self) -> tuple[argparse.ArgumentParser, LoggerConfig]:
        """Return the argument parser and default logger configuration."""

        raise NotImplementedError

    def prepare_arguments(
        self,
        parser: argparse.ArgumentParser,
        args: argparse.Namespace,
        argv: Sequence[str] | None,
    ) -> argparse.Namespace:
        """Hook executed after parsing command line arguments."""

        return args

    def handle_pre_run(
        self, parser: argparse.ArgumentParser, args: argparse.Namespace
    ) -> int | None:
        """Return an exit code to abort execution before running the pipeline."""

        return None

    def get_program_name(self) -> str:
        """Return the identifier used when creating log files."""

        return Path(__file__).with_suffix("").name

    def get_logging_date(self, args: argparse.Namespace) -> str | None:
        """Return the date token forwarded to :func:`setup_cli_logging`."""

        return getattr(args, "date", None)

    def get_logger(self) -> Logger:
        """Return the logger instance used for pipeline level events."""

        return default_logger

    def get_config_mapping(self) -> Mapping[str, str]:
        """Return CLI-to-config attribute mappings forwarded to the loader."""

        return {}

    def get_base_parser(
        self, parser: argparse.ArgumentParser
    ) -> argparse.ArgumentParser | None:
        """Return the root parser forwarded to :func:`run_cli_command`."""

        del parser
        return None

    def run_pipeline(self, cfg: Any, args: argparse.Namespace) -> int:
        """Execute the pipeline logic for the concrete implementation."""

        raise NotImplementedError

    def on_logging_ready(self, logging_ctx: CLILoggingContext) -> None:
        """Hook executed after log handlers have been initialised."""

    def execute(
        self,
        args: argparse.Namespace,
        parser: argparse.ArgumentParser,
        logging_ctx: CLILoggingContext,
    ) -> int:
        """Execute the pipeline using :func:`run_cli_command`."""

        from ..cli_utils import run_cli_command

        return run_cli_command(
            args=args,
            parser=parser,
            base_parser=self.get_base_parser(parser),
            log_cfg=logging_ctx.log_cfg,
            mapping=dict(self.get_config_mapping()),
            run=self.run_pipeline,
            logger=self.get_logger(),
        )

    def after_run(self, log_cfg: LoggerConfig, exit_code: int) -> int:
        """Finalisation hook executed once the CLI run finishes."""

        configure_logger(log_cfg)
        return exit_code

    def main(self, argv: Sequence[str] | None = None) -> int:
        """Entry point mirroring the previous function-based CLI contract."""

        parser, log_cfg = self.build_parser()
        args = parser.parse_args(argv)
        args = self.prepare_arguments(parser, args, argv)
        exit_code = self.handle_pre_run(parser, args)
        if exit_code is not None:
            return exit_code

        program_name = self.get_program_name()
        date_token = self.get_logging_date(args)

        with setup_cli_logging(program_name, log_cfg, date_token) as logging_ctx:
            self.on_logging_ready(logging_ctx)
            exit_code = self.execute(args, parser, logging_ctx)

        return self.after_run(log_cfg, exit_code)


__all__ = ["PipelineCLIBase"]
