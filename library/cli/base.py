"""Base classes for command-line pipeline entry points."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

from ..common.log import logger as default_logger
from .logging import CLILoggingContext, setup_cli_logging
from .parser import Logger, LoggerConfig, configure_logger
from .run_context import compute_generated_at


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

        module = sys.modules.get(self.__class__.__module__)
        module_file = getattr(module, "__file__", None) if module is not None else None
        if module_file:
            return Path(module_file).with_suffix("").name
        return self.__class__.__name__.lower()

    def get_logging_date(self, args: argparse.Namespace) -> str | None:
        """Return the date token forwarded to :func:`setup_cli_logging`."""

        return getattr(args, "date", None)

    def resolve_generated_at(
        self,
        args: argparse.Namespace,
        argv: Sequence[str] | None,
        log_cfg: LoggerConfig,
        *,
        date_token: str | None,
    ) -> str:
        """Return a deterministic ``generated_at`` value for metadata outputs."""

        invocation = getattr(args, "invocation", None)
        seed_parts: list[str] = []
        if isinstance(invocation, Sequence) and invocation:
            seed_parts.extend(str(part) for part in invocation)
        else:
            program = self.get_program_name()
            seed_parts.append(program)
            if argv is not None:
                seed_parts.extend(str(part) for part in argv)

        return compute_generated_at(
            date_token=date_token,
            run_id=log_cfg.run_id,
            seed_parts=seed_parts,
        )

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

    def get_run_cli_command(self) -> Callable[..., int]:
        """Return the :func:`run_cli_command` implementation to use."""

        module = sys.modules.get(self.__class__.__module__)
        if module is not None:
            candidate = getattr(module, "run_cli_command", None)
            if callable(candidate):
                return candidate
        from ..cli_utils import run_cli_command as _run_cli_command

        return _run_cli_command

    def execute(
        self,
        args: argparse.Namespace,
        parser: argparse.ArgumentParser,
        logging_ctx: CLILoggingContext,
    ) -> int:
        """Execute the pipeline using the resolved :func:`run_cli_command`."""

        run_cli_command = self.get_run_cli_command()

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

        explicit_run_id = getattr(args, "run_id", None)
        if isinstance(explicit_run_id, str):
            explicit_run_id = explicit_run_id.strip() or None
        if explicit_run_id in (argparse.SUPPRESS,):
            explicit_run_id = None
        args.run_id = explicit_run_id
        if explicit_run_id is not None:
            log_cfg.run_id = explicit_run_id

        if not hasattr(args, "invocation"):
            from ..cli_utils import resolve_invocation

            invocation = resolve_invocation(parser.prog, argv)
            args.invocation = invocation

        run_id_value = getattr(args, "run_id", None)
        if not run_id_value and not log_cfg.run_id:
            from ..cli_utils import _canonical_run_descriptor, uuid

            descriptor = _canonical_run_descriptor(args, parser)
            if descriptor:
                run_id_value = uuid.uuid5(uuid.NAMESPACE_URL, descriptor).hex
                args.run_id = run_id_value
                log_cfg.run_id = run_id_value

        exit_code = self.handle_pre_run(parser, args)
        if exit_code is not None:
            return exit_code

        program_name = self.get_program_name()
        date_token = self.get_logging_date(args)
        log_cfg.generated_at = self.resolve_generated_at(
            args,
            argv,
            log_cfg,
            date_token=date_token,
        )

        with setup_cli_logging(program_name, log_cfg, date_token) as logging_ctx:
            self.on_logging_ready(logging_ctx)
            exit_code = self.execute(args, parser, logging_ctx)

        return self.after_run(log_cfg, exit_code)


__all__ = ["PipelineCLIBase", "compute_generated_at"]
