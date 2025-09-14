"""Command line interface for retrieving ChEMBL activity data."""

from __future__ import annotations

import sys

# ruff: noqa: E402
from pathlib import Path

if __package__ is None:  # running as a script
    sys.path.append(str(Path(__file__).resolve().parents[1]))

import argparse
from collections.abc import Iterable, Sequence
from itertools import islice

import requests
from pandera.errors import SchemaErrors

from library import chembl_library as cl
from library import io
from library.chembl_client import ChemblClient
from library.cli import (
    LoggerConfig,
    apply_config_overrides,
    configure_logger,
)
from library.cli import (
    build_parser as base_parser,
)
from library.config import (
    Config,
    _serialize_paths,
    ensure_dirs,
    print_config,
)
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from schemas import ActivitiesSchema, normalize_activities


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute activity retrieval from the ChEMBL API.

    The resulting CSV places columns defined in :data:`ActivitiesSchema`
    first, preserving their declared order. Any additional fields appear
    afterwards sorted alphabetically.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args : argparse.Namespace
        Parsed command-line arguments. ``args.limit`` constrains the number of
        identifiers processed, while ``args.dry_run`` skips network calls and
        file generation.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    limit = cfg.activity.limit
    if limit is not None and limit < 0:
        logger.error("activity.limit must be non-negative")
        return 1

    if cfg.activity.dry_run:
        expected = limit if limit is not None else 0
        logger.info("dry_run", limit=expected)
        return 0

    # Configure HTTP session with the supplied User-Agent and retry policy
    with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
        try:
            ids_iter = io.read_ids(
                args.input_csv, column=cfg.activity.column, cfg=cfg.io
            )
        except (FileNotFoundError, ValueError) as exc:
            logger.error("%s", exc)
            return 1

        # Apply the ``limit`` without materialising the entire iterator first.
        # ``itertools.islice`` allows lazy slicing; converting to ``list`` enables
        # length calculation for logging purposes.
        ids: Iterable[str] = ids_iter
        if limit is not None:
            limited = list(islice(ids_iter, limit))
            ids = limited
            logger.info("process_limit", limit=len(limited))

        try:
            df = cl.get_activities(
                ids,
                cfg=cfg.api,
                client=client,
                chunk_size=cfg.activity.chunk_size,
                timeout=cfg.activity.timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error("failed to retrieve activities: %s", exc)
            return 1
        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        df = normalize_activities(df)
        # Determine final column order: schema-defined columns first in their
        # declared sequence, followed by any additional columns sorted
        # alphabetically to provide deterministic output.
        schema_cols = list(ActivitiesSchema.columns)
        head = [c for c in schema_cols if c in df.columns]
        tail = sorted(c for c in df.columns if c not in schema_cols)
        col_order = head + tail
        rows_total = len(df)
        exit_code = 0
        required_cols = {
            name for name, col in ActivitiesSchema.columns.items() if col.required
        }
        optional_cols = set(ActivitiesSchema.columns) - required_cols
        missing_required = required_cols - set(df.columns)
        missing_optional = optional_cols - set(df.columns)
        if not missing_required:
            if missing_optional:
                logger.warning(
                    "DataFrame is missing optional columns: %s", missing_optional
                )
            try:
                df = ActivitiesSchema.validate(df, lazy=True)
            except SchemaErrors as exc:
                failure_path = output.with_name(f"{output.stem}_failure_cases.csv")
                errors = SidecarErrors()
                for row in exc.failure_cases.to_dict("records"):
                    errors.add_error(row)
                errors.save(failure_path)
                logger.error(
                    "validation failed; wrote %d failure cases to %s",
                    len(exc.failure_cases),
                    failure_path,
                )
                df = getattr(exc, "validated_data", df)
                exit_code = 1
        else:
            logger.warning(
                "Skipping validation due to missing required columns: %s",
                missing_required,
            )
        rows_kept = len(df)
        rows_dropped = rows_total - rows_kept
        try:
            key_cols = [c for c in ["activity_id"] if c in df.columns]
            csv_path = io.write_csv(
                df,
                output,
                cfg=cfg,
                key_cols=key_cols or None,
                col_order=col_order,
            )
            logger.info("write_done", rows=rows_kept, path=str(csv_path))
        except OSError as exc:
            logger.error("failed to write output CSV: %s", exc)
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
            schema="ActivitiesSchema",
        )

        try:
            analyze_table_quality(df, table_name=str(output.with_suffix("")))
        except ValueError as exc:
            logger.error("failed to generate quality report: %s", exc)
            return 1
        return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser, log_cfg = base_parser(
        "ChEMBL activity data utilities", column="activity_id", chunk_size=5
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
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Read input and exit without contacting the API or writing files",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit <= 0:
        # Reject non-positive limits early to provide clear CLI feedback.
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
                "timeout": "activity.timeout",
                "column": "activity.column",
                "chunk_size": "activity.chunk_size",
                "limit": "activity.limit",
                "dry_run": "activity.dry_run",
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
        logger.error("%s", exc)
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("failed to set up directories: %s", exc)
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
