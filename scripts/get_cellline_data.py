from __future__ import annotations

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from . import _bootstrap as _bootstrap_module
elif __package__ in {None, ""}:
    import _bootstrap as _bootstrap_module  # pragma: no cover - CLI fallback
else:  # pragma: no cover - executed when imported as a package module
    from . import _bootstrap as _bootstrap_module

bootstrap_cli = _bootstrap_module.bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli
del _bootstrap_module

import pandas as pd

from library import cli, io
from library.cli import LoggerConfig, configure_logger
from library.cli import build_parser as base_parser
from library.cli.logging import setup_cli_logging
from library.cli_utils import run_cli_command
from library.common.log import logger
from library.config import Config
from library.pipelines.cellline import (
    CellLinePipelineOptions,
    run_cellline_pipeline,
)
from library.utils.data_correlation import generate_correlation_report
from library.utils.qc_report import generate_qc_report

DEFAULT_INPUT_NAME = "cellline.csv"
DEFAULT_OUTPUT_STEM = "cellline"
MODE_CHOICES: tuple[str, ...] = ("chembl",)


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line parser for the cell line pipeline."""

    parser, log_cfg = base_parser(
        "ChEMBL cell line data utilities",
        column="cell_chembl_id",
        chunk_size=20,
        size_option="--batch-size",
        size_dest="batch_size",
    )
    parser.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
    parser.add_argument(
        "--mode",
        choices=MODE_CHOICES,
        default="chembl",
        help="Data source to query (default: chembl)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=None,
        help="Override request timeout in seconds",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process (0 skips processing)",
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    parser.set_defaults(func=run)
    return parser, log_cfg


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the cell line pipeline using ``cfg`` and parsed ``args``."""

    mode = str(getattr(args, "mode", "chembl")).lower()
    if mode not in MODE_CHOICES:
        logger.error("cellline_unsupported_mode", mode=mode)
        return 1

    output_path = Path(
        args.output_csv
        or io.default_output_path(
            args.input_csv,
            cfg.io,
            date=getattr(args, "date", None),
        )
    )
    args.output_csv = output_path

    if (
        getattr(args, "skip_existing", False)
        and output_path.exists()
        and not getattr(args, "force", False)
    ):
        logger.info("pipeline_skip_existing", output=str(output_path))
        return 0

    batch_size = getattr(args, "batch_size", cfg.cellline.batch_size)
    timeout = args.timeout if args.timeout is not None else None

    options = CellLinePipelineOptions(
        input_csv=Path(args.input_csv),
        output_csv=output_path,
        column=str(args.column),
        batch_size=int(batch_size),
        limit=args.limit,
        offset=int(args.offset),
        timeout=timeout,
        mode=mode,
    )

    logger.info(
        "cellline_pipeline_start",
        input=str(args.input_csv),
        output=str(output_path),
        limit=args.limit,
        offset=args.offset,
        timeout=timeout if timeout is not None else cfg.cellline.timeout,
        batch_size=options.batch_size,
    )

    result = run_cellline_pipeline(cfg, options)

    exit_code = result.exit_code
    artifacts: io.StandardOutputArtifacts | None = None
    output_path = result.output_path

    if exit_code == 0 and result.written:
        try:
            dataset_frame = pd.read_csv(
                output_path,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
            )
        except (
            OSError,
            ValueError,
            pd.errors.EmptyDataError,
            pd.errors.ParserError,
        ) as exc:
            logger.error(
                "cellline_output_load_failed",
                error=str(exc),
                path=str(output_path),
            )
            exit_code = 1
        else:
            table_name, date_tag = io.derive_output_labels(
                output_path,
                default_table=DEFAULT_OUTPUT_STEM,
            )

            doc_quality_cfg = getattr(getattr(cfg, "system", None), "doc_quality", None)
            include_columns = (
                getattr(doc_quality_cfg, "include_columns", None)
                if doc_quality_cfg is not None
                else None
            )
            exclude_columns = (
                getattr(doc_quality_cfg, "exclude_columns", None)
                if doc_quality_cfg is not None
                else None
            )
            sample_rows = (
                getattr(doc_quality_cfg, "sample_rows", None)
                if doc_quality_cfg is not None
                else None
            )
            profiler = (
                getattr(doc_quality_cfg, "profiler", None)
                if doc_quality_cfg is not None
                else None
            )

            try:
                correlation_report = generate_correlation_report(
                    dataset_frame,
                    table_name=table_name,
                    include_columns=include_columns,
                    exclude_columns=exclude_columns,
                    sample_rows=sample_rows,
                    profiler=profiler,
                )
            except Exception as exc:  # pragma: no cover - defensive logging
                logger.warning(
                    "cellline_correlation_generation_failed",
                    error=str(exc),
                    path=str(output_path),
                )
                correlation_report = pd.DataFrame()

            try:
                quality_report = generate_qc_report(
                    dataset_frame,
                    table_name=table_name,
                    include_columns=include_columns,
                    exclude_columns=exclude_columns,
                    sample_rows=sample_rows,
                    profiler=profiler,
                )
            except Exception as exc:  # pragma: no cover - defensive logging
                logger.warning(
                    "cellline_quality_generation_failed",
                    error=str(exc),
                    path=str(output_path),
                )
                quality_report = pd.DataFrame()

            try:
                artifacts = io.save_standard_outputs(
                    dataset_frame,
                    correlation_report,
                    quality_report,
                    table_name=table_name,
                    date_tag=date_tag,
                    key_columns=["cell_chembl_id"],
                    output_path=output_path,
                )
            except (OSError, ValueError) as exc:
                logger.error(
                    "cellline_standard_outputs_failed",
                    error=str(exc),
                    path=str(output_path),
                )
                exit_code = 1
            else:
                logger.info(
                    "cellline_standard_outputs_written",
                    dataset=str(artifacts.dataset),
                    quality_report=str(artifacts.quality_report),
                    correlation_report=str(artifacts.correlation_report),
                )
                args.output_csv = artifacts.dataset
                args._cellline_artifacts = artifacts
                output_path = artifacts.dataset

                try:
                    metadata_path = io.save_metadata(
                        table_name=table_name,
                        date_tag=date_tag,
                        args=args,
                        qc_summary={"rows": int(len(dataset_frame))},
                        output_dir=artifacts.dataset.parent,
                        artifacts=[
                            artifacts.dataset,
                            artifacts.quality_report,
                            artifacts.correlation_report,
                        ],
                    )
                except (OSError, ValueError) as exc:
                    logger.error(
                        "cellline_metadata_persist_failed",
                        error=str(exc),
                        output_dir=str(artifacts.dataset.parent),
                    )
                    exit_code = 1
                else:
                    logger.info(
                        "cellline_metadata_written",
                        metadata=str(metadata_path),
                    )

    preserve_intermediate = any(
        bool(getattr(args, flag, False))
        for flag in ("debug", "keep_intermediate", "emit_legacy_artifacts")
    )

    if result.failure_path and result.failures:
        logger.error(
            "cellline_validation_failed",
            failures=result.failures,
            path=str(result.failure_path),
        )
        if not preserve_intermediate:
            result.failure_path.unlink(missing_ok=True)
    elif result.failure_path and not preserve_intermediate:
        result.failure_path.unlink(missing_ok=True)

    if result.missing_ids:
        sample = list(result.missing_ids[:5])
        logger.warning(
            "cellline_missing_identifiers_summary",
            total=len(result.missing_ids),
            sample=sample,
        )

    status = "OK" if exit_code == 0 else "ERROR"
    logger.info(
        "cellline_pipeline_summary",
        records=result.records,
        duration=result.duration,
        output=str(output_path),
        exit_code=exit_code,
        status=status,
    )
    return exit_code


def main(argv: Sequence[str] | None = None) -> int:
    """Program entry point for the cell line CLI."""

    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    cli.prepare_io_paths(
        args,
        input_default=DEFAULT_INPUT_NAME,
        output_stem=DEFAULT_OUTPUT_STEM,
    )
    if args.limit == 0:
        logger.info("pipeline_skip_limit", limit=args.limit)
        return 0
    if args.limit is not None and args.limit < 0:
        parser.error("--limit must be zero or a positive integer")
    if args.offset < 0:
        parser.error("--offset must be zero or a positive integer")

    with setup_cli_logging(
        Path(__file__).with_suffix("").name,
        log_cfg,
        getattr(args, "date", None),
    ) as logging_ctx:
        exit_code = run_cli_command(
            args=args,
            parser=parser,
            log_cfg=logging_ctx.log_cfg,
            mapping={
                "timeout": "cellline.timeout",
                "column": "cellline.column",
                "batch_size": "cellline.batch_size",
                "limit": "cellline.limit",
                "offset": "cellline.offset",
            },
            run=run,
            logger=logger,
        )
    configure_logger(log_cfg)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
