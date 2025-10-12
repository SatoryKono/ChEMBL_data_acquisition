"""Standalone CLI tools used during development and QA."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Mapping, Sequence

from library.cli import Logger, LoggerConfig
from library.cli_utils import ensure_run_id, resolve_invocation, run_cli_command
from library.config import Config

PrepareCallback = Callable[
    [argparse.ArgumentParser, argparse.Namespace, Sequence[str] | None],
    argparse.Namespace | None,
]


def run_cli_tool(
    *,
    build_parser: Callable[[], tuple[argparse.ArgumentParser, LoggerConfig]],
    run: Callable[[Config, argparse.Namespace], int],
    argv: Sequence[str] | None = None,
    mapping: Mapping[str, str] | None = None,
    logger: Logger | None = None,
    base_parser: argparse.ArgumentParser | None = None,
    prepare: PrepareCallback | None = None,
) -> int:
    """Execute a CLI tool using :func:`library.cli_utils.run_cli_command`.

    Parameters
    ----------
    build_parser:
        Callable returning the command specific parser and its logging
        configuration.
    run:
        Callable implementing the command behaviour. It receives the resolved
        :class:`~library.config.Config` and the parsed arguments.
    argv:
        Optional sequence of command line arguments forwarded to the parser.
        ``None`` falls back to :data:`sys.argv[1:]` mirroring
        :mod:`argparse` semantics.
    mapping:
        Mapping of CLI argument names to configuration paths passed to
        :func:`library.cli.parser.apply_config_overrides`.
    logger:
        Optional logger instance used for structured logging. When omitted the
        default logger configured by :func:`library.cli.configure_logger` is
        used.
    base_parser:
        Optional parser providing default values for shared arguments. This is
        forwarded to :func:`library.cli_utils.run_cli_command`.
    prepare:
        Optional callback invoked after parsing. It can adjust ``args`` prior
        to delegating to :func:`library.cli_utils.run_cli_command`.

    Returns
    -------
    int
        Process exit code returned by :func:`library.cli_utils.run_cli_command`.
    """

    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)

    if prepare is not None:
        prepared = prepare(parser, args, argv)
        if prepared is not None:
            args = prepared

    if not hasattr(args, "invocation"):
        args.invocation = resolve_invocation(parser.prog, argv)

    ensure_run_id(args, parser, log_cfg)

    return run_cli_command(
        args=args,
        parser=parser,
        base_parser=base_parser,
        log_cfg=log_cfg,
        mapping=dict(mapping or {}),
        run=run,
        logger=logger,
    )


__all__ = ["run_cli_tool"]
