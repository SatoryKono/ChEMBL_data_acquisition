"""Command line interface for retrieving ChEMBL assay data.

The script provides a reusable :func:`main` entry point together with helper
functions that unit tests import directly. Functions return explicit exit codes
instead of terminating the interpreter to make orchestration easier.
"""

from __future__ import annotations

import sys
from pathlib import Path

import argparse
from collections.abc import Iterable, Sequence
from functools import partial
from itertools import islice

import pandas as pd
import requests

# Ensure repository package imports work when the script is executed directly.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from library import assay_postprocessing as ap
from library import chembl_library as cl
from library import cli
from library import io
from library.clients import ChemblClient
from library.cli import (
    LoggerConfig,
    configure_logger,
)
from library.cli import build_parser as base_parser
from library.cli_utils import PipelineError, run_pipeline
from library.config import Config, _serialize_paths, ensure_dirs, print_config
from library.log import logger
from library.pipeline_metadata import add_pipeline_metadata
from library.table_quality import analyze_table_quality
from library.validation import validate_assays
from schemas import AssaysSchema, normalize_assays

__all__ = ["ap", "main"]


DEFAULT_INPUT_NAME = "assay.csv"
DEFAULT_OUTPUT_STEM = "assays"


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute assay retrieval from the ChEMBL API.

    The output CSV arranges columns so that fields defined in
    :class:`~schemas.assays.AssaysSchema` appear first. Any additional columns
    are appended alphabetically.

    Parameters
    ----------
    cfg : Config
        Application configuration providing API credentials, retry strategy and
        CSV export settings.
    args : argparse.Namespace
        Parsed command-line arguments accepted by the ``assay`` CLI command.

    Returns
    -------
    int
        ``0`` on success, non-zero when validation or I/O failures occur.
        Upstream API errors are logged and converted into a failure code by
        :func:`library.cli_utils.run_pipeline`.
    """
    limit = cfg.assay.limit
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="assay.limit", limit=limit)
        return 1

    try:
        ids_iter = io.read_ids(args.input_csv, column=cfg.assay.column, cfg=cfg.io)
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "read_fail",
            error=str(exc),
            path=str(args.input_csv),
        )
        return 1

    offset = getattr(args, "offset", 0)
    if offset:
        ids_iter = islice(ids_iter, offset, None)
        logger.info("process_offset", offset=offset)

    ids_source: Iterable[str]
    if limit is not None:
        limited_ids = list(islice(ids_iter, limit))
        ids_source = limited_ids
        logger.info("process_limit", limit=len(limited_ids))
    else:
        ids_source = ids_iter

    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    failure_path = Path(output).with_name(f"{Path(output).stem}_failure_cases.csv")

    def fetcher() -> Iterable[pd.DataFrame]:
        with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
            try:
                df = cl.get_assays(
                    ids_source,
                    cfg=cfg.api,
                    client=client,
                    chunk_size=cfg.assay.batch_size,
                    timeout=cfg.assay.timeout,
                )
            except (requests.RequestException, ValueError) as exc:
                logger.error(
                    "assay_fetch_failed",
                    extra={"msg": str(exc)},
                    error=str(exc),
                    batch_size=cfg.assay.batch_size,
                    timeout=cfg.assay.timeout,
                )
                raise PipelineError(str(exc)) from exc
            yield df

    metadata_hooks = [
        ap.postprocess_assays,
        normalize_assays,
        add_pipeline_metadata,
    ]

    validators = [partial(validate_assays, return_result=True)]

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: Sequence[str],
        key_cols: Sequence[str],
    ) -> Path:
        frames = [chunk for chunk in chunks]
        if frames:
            df = pd.concat(frames, ignore_index=True)
        else:
            df = pd.DataFrame(columns=col_order)
        resolved_keys = list(key_cols) or None
        return io.write_csv(
            df,
            destination,
            cfg=cfg,
            key_cols=resolved_keys,
            col_order=col_order,
        )

    table_quality = partial(
        analyze_table_quality,
        table_name=str(Path(output).with_suffix("")),
    )

    return run_pipeline(
        fetcher=fetcher,
        schema=AssaysSchema,
        schema_name="AssaysSchema",
        validators=validators,
        metadata_hooks=metadata_hooks,
        writer=writer,
        output_path=output,
        failure_path=failure_path,
        command=" ".join(sys.argv),
        config_snapshot=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(args.input_csv)},
        key_columns=["assay_chembl_id"],
        table_quality=table_quality,
        logger=logger,
    )


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        A tuple containing the fully configured parser and default logging
        configuration for the pipeline run.
    """
    parser, log_cfg = base_parser(
        "ChEMBL assay data utilities",
        column="assay_chembl_id",
        chunk_size=10,
        size_option="--batch-size",
        size_dest="batch_size",
    )
    parser.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the assay pipeline with optional argument overrides.

    Parameters
    ----------
    argv : Sequence[str] | None, optional
        Command-line arguments to parse. When ``None`` the process arguments
        from :data:`sys.argv` are used.

    Returns
    -------
    int
        ``0`` on success, non-zero otherwise. Failures are logged with context
        describing the failing section.

    Raises
    ------
    SystemExit
        Raised when invalid command-line options are supplied.
    """
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    cli.prepare_io_paths(
        args,
        input_default=DEFAULT_INPUT_NAME,
        output_stem=DEFAULT_OUTPUT_STEM,
    )
    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be a positive integer")
    if args.offset < 0:
        parser.error("--offset must be zero or a positive integer")
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = cli.apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "timeout": "assay.timeout",
                "column": "assay.column",
                "batch_size": "assay.batch_size",
                "limit": "assay.limit",
            },
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        output_path = Path(
            args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        )
        args.output_csv = output_path
        if args.skip_existing and output_path.exists() and not args.force:
            logger.info("pipeline_skip_existing", output=str(output_path))
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        logger = configure_logger(log_cfg)
    except (ValueError, TypeError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error(
            "directory_setup_failed",
            error=str(exc),
            output=str(args.output_csv),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    exit_code: int = args.func(cfg, args)
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
