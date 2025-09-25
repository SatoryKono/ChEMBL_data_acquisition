"""Command line interface for retrieving ChEMBL assay data."""

from __future__ import annotations

import sys

# ruff: noqa: E402
from pathlib import Path

if __package__ is None:  # running as a script
    sys.path.append(str(Path(__file__).resolve().parents[1]))

import argparse
from collections.abc import Sequence
from itertools import islice

import requests
from pandera.errors import SchemaErrors

from library import assay_postprocessing as ap
from library import chembl_library as cl
from library import io
from library.chembl_client import ChemblClient
from library.cli import (
    LoggerConfig,
    apply_config_overrides,
    configure_logger,
)
from library.cli import build_parser as base_parser
from library.config import (
    Config,
    _serialize_paths,
    ensure_dirs,
    print_config,
)
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.pipeline_metadata import add_pipeline_metadata
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from library.validation import validate_assays
from schemas import AssaysSchema, normalize_assays

__all__ = ["ap", "main"]


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute assay retrieval from the ChEMBL API.

    The output CSV arranges columns so that fields defined in
    :class:`~schemas.assays.AssaysSchema` appear first.  Any additional columns
    are appended alphabetically.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    limit = cfg.assay.limit
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="assay.limit", limit=limit)
        return 1

    # Prepare HTTP session for ChEMBL requests
    with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
        try:
            ids_iter = io.read_ids(args.input_csv, column=cfg.assay.column, cfg=cfg.io)
        except (FileNotFoundError, ValueError) as exc:
            logger.error(
                "read_fail",
                error=str(exc),
                path=str(args.input_csv),
            )
            return 1

        ids = ids_iter
        if limit is not None:
            limited_ids = list(islice(ids_iter, limit))
            ids = limited_ids
            logger.info("process_limit", limit=len(limited_ids))

        try:
            df = cl.get_assays(
                ids,
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
            return 1
        df = ap.postprocess_assays(df)
        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        df = normalize_assays(df)
        df = add_pipeline_metadata(df)
        rows_total = len(df)
        exit_code = 0
        required_cols = {
            name for name, col in AssaysSchema.columns.items() if col.required
        }
        optional_cols = set(AssaysSchema.columns) - required_cols
        missing_required = required_cols - set(df.columns)
        missing_optional = optional_cols - set(df.columns)
        if not missing_required:
            if missing_optional:
                logger.warning(
                    "optional_columns_missing",
                    columns=sorted(missing_optional),
                )
            try:
                validation_result = validate_assays(df, return_result=True)
            except SchemaErrors as exc:
                failure_path = Path(output).with_name(
                    f"{Path(output).stem}_failure_cases.csv"
                )
                errors = SidecarErrors()
                for row in exc.failure_cases.to_dict("records"):
                    errors.add_error(row)
                errors.save(failure_path)
                logger.error(
                    "validation_failed",
                    failures=len(exc.failure_cases),
                    path=str(failure_path),
                )
                df = getattr(exc, "validated_data", df)
                exit_code = 1
            else:
                df = validation_result.data
                if not validation_result.failure_cases.empty:
                    failure_path = Path(output).with_name(
                        f"{Path(output).stem}_failure_cases.csv"
                    )
                    errors = SidecarErrors()
                    for row in validation_result.failure_cases.to_dict("records"):
                        errors.add_error(row)
                    errors.save(failure_path)
                    logger.error(
                        "validation_failed",
                        failures=len(validation_result.failure_cases),
                        path=str(failure_path),
                    )
                    exit_code = 1
        else:
            logger.warning(
                "validation_skipped",
                missing_columns=sorted(missing_required),
            )
        rows_kept = len(df)
        rows_dropped = rows_total - rows_kept
        # Arrange columns so schema-defined fields appear first while any
        # additional columns are sorted alphabetically and placed at the end.
        schema_cols: list[str] = list(AssaysSchema.columns)
        head = [c for c in schema_cols if c in df.columns]
        tail = sorted(c for c in df.columns if c not in schema_cols)
        col_order = head + tail
        try:
            key_cols = [c for c in ["assay_chembl_id"] if c in df.columns]
            csv_path = io.write_csv(
                df,
                output,
                cfg=cfg,
                key_cols=key_cols or None,
                col_order=col_order,
            )
            logger.info("write_done", rows=rows_kept, path=str(csv_path))
        except OSError as exc:
            logger.error(
                "write_fail",
                error=str(exc),
                path=str(output),
            )
            return 1

        stats: Stats = {
            "rows_total": rows_total,
            "rows_kept": rows_kept,
            "rows_dropped": rows_dropped,
            "output_sha256": file_sha256(csv_path),
        }
        write_meta_yaml(
            csv_path=csv_path,
            command=" ".join(sys.argv),
            config_subset=_serialize_paths(cfg.to_dict()),
            inputs={"input_csv": str(args.input_csv)},
            stats=stats,
            schema="AssaysSchema",
        )

        try:
            analyze_table_quality(df, table_name=str(output.with_suffix("")))
        except ValueError as exc:
            logger.error(
                "quality_report_failed",
                error=str(exc),
                path=str(output),
            )
            return 1
        return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser, log_cfg = base_parser(
        "ChEMBL assay data utilities",
        column="assay_chembl_id",
        chunk_size=10,
        size_option="--batch-size",
        size_dest="batch_size",
    )
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
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be a positive integer")
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = apply_config_overrides(
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
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
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
