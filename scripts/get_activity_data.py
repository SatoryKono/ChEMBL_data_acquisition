"""Command line interface for retrieving ChEMBL activity data."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterator, Sequence
from itertools import islice
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd


def _ensure_project_root() -> None:
    """Ensure the repository root is on ``sys.path`` when executed as a script."""

    script_path = Path(__file__).resolve()
    project_root = script_path.parents[1]
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)


if __package__ in {None, ""}:
    _ensure_project_root()

import requests
from pandera.errors import SchemaErrors

from library import chembl_library as cl
from library import cli
from library import io
from library.csv_utils import write_csv_chunks_deterministic
from library.clients import ChemblClient, _chunked
from library.cli import (
    LoggerConfig,
    configure_logger,
)
from library.cli import (
    build_parser as base_parser,
)
from library.config import (
    ActivityActionTypeCfg,
    ActivityBoundsCfg,
    ActivityPropertiesCfg,
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
from library.validation import validate_activities
from library.processing.activity import (
    apply_activity_annotations,
    compute_activity_bounds,
)
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
        logger.error("invalid_limit", section="activity.limit", limit=limit)
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

        # Apply the ``limit`` without materialising the entire iterator first.
        ids = islice(ids_iter, limit) if limit is not None else ids_iter

        enrichment_cfg = cfg.activity_enrichment
        extra_columns: list[str] = []
        action_cfg = enrichment_cfg.action_type
        if action_cfg.enabled or action_cfg.log_missing or action_cfg.log_distribution:
            extra_columns.append(action_cfg.column)
        extra_kwargs = {"extra_columns": extra_columns} if extra_columns else {}

        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)

        processed_ids = 0

        def _iter_ids() -> Iterator[str]:
            nonlocal processed_ids
            for identifier in ids:
                processed_ids += 1
                yield identifier

        id_chunks = _chunked(_iter_ids(), cfg.activity.batch_size)

        required_cols = {
            name for name, col in ActivitiesSchema.columns.items() if col.required
        }
        optional_cols = set(ActivitiesSchema.columns) - required_cols
        present_columns: set[str] = set()
        total_failures = 0
        rows_total = 0
        rows_kept = 0
        rows_dropped = 0
        exit_code = 0
        failure_path = Path(output).with_name(f"{Path(output).stem}_failure_cases.csv")
        errors = SidecarErrors()

        chunk_paths: list[Path] = []
        all_columns: set[str] = set()

        validation_enabled = True
        missing_required_columns: set[str] = set()

        csv_path: Path | None = None
        with TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            for chunk_index, chunk_ids in enumerate(id_chunks):
                try:
                    chunk_df = cl.get_activities(
                        chunk_ids,
                        cfg=cfg.api,
                        client=client,
                        chunk_size=cfg.activity.batch_size,
                        timeout=cfg.activity.timeout,
                        **extra_kwargs,
                    )
                except (requests.RequestException, ValueError) as exc:
                    logger.error(
                        "activity_fetch_failed",
                        extra={"msg": str(exc)},
                        error=str(exc),
                        batch_size=cfg.activity.batch_size,
                        timeout=cfg.activity.timeout,
                    )
                    return 1

                chunk_df = normalize_activities(chunk_df)
                chunk_df = add_pipeline_metadata(chunk_df)
                chunk_df = compute_activity_bounds(chunk_df, cfg.activity_bounds)
                chunk_df = apply_activity_annotations(
                    chunk_df,
                    action_cfg=enrichment_cfg.action_type,
                    properties_cfg=enrichment_cfg.activity_properties,
                )

                chunk_columns = set(chunk_df.columns)
                present_columns.update(chunk_columns)
                missing_chunk_required = required_cols - chunk_columns
                if missing_chunk_required:
                    missing_required_columns.update(missing_chunk_required)
                    validation_enabled = False

                rows_total += len(chunk_df)

                validated_chunk = chunk_df
                if validation_enabled and not chunk_df.empty:
                    try:
                        validation_result = validate_activities(
                            chunk_df, return_result=True
                        )
                    except SchemaErrors as exc:
                        for row in exc.failure_cases.to_dict("records"):
                            errors.add_error(row)
                        total_failures += len(exc.failure_cases)
                        logger.error(
                            "validation_failed",
                            failures=len(exc.failure_cases),
                            path=str(failure_path),
                        )
                        validated_chunk = getattr(exc, "validated_data", chunk_df)
                        exit_code = 1
                    else:
                        validated_chunk = validation_result.data
                        if not validation_result.failure_cases.empty:
                            cases = validation_result.failure_cases.to_dict("records")
                            for row in cases:
                                errors.add_error(row)
                            total_failures += len(cases)
                            logger.error(
                                "validation_failed",
                                failures=len(cases),
                                path=str(failure_path),
                            )
                            exit_code = 1

                rows_kept += len(validated_chunk)
                rows_dropped += len(chunk_df) - len(validated_chunk)
                present_columns.update(validated_chunk.columns)
                all_columns.update(validated_chunk.columns)

                chunk_path = Path(tmpdir) / f"chunk_{chunk_index}.pkl"
                validated_chunk.to_pickle(chunk_path)
                chunk_paths.append(chunk_path)

            if not chunk_paths:
                empty_chunk = apply_activity_annotations(
                    compute_activity_bounds(
                        add_pipeline_metadata(pd.DataFrame()), cfg.activity_bounds
                    ),
                    action_cfg=enrichment_cfg.action_type,
                    properties_cfg=enrichment_cfg.activity_properties,
                )
                present_columns.update(empty_chunk.columns)
                all_columns.update(empty_chunk.columns)
                chunk_path = Path(tmpdir) / "chunk_0.pkl"
                empty_chunk.to_pickle(chunk_path)
                chunk_paths.append(chunk_path)

            schema_cols = list(ActivitiesSchema.columns)
            head = [c for c in schema_cols if c in all_columns]
            tail = sorted(c for c in all_columns if c not in schema_cols)
            col_order = head + tail

            def _iter_validated_chunks() -> Iterator[pd.DataFrame]:
                for path in chunk_paths:
                    df_chunk = pd.read_pickle(path)
                    if col_order:
                        df_chunk = df_chunk.reindex(columns=col_order)
                    yield df_chunk

            try:
                key_cols = [c for c in ["activity_id"] if c in col_order]
                sort_columns = key_cols or sorted(col_order)
                csv_path = write_csv_chunks_deterministic(
                    _iter_validated_chunks(),
                    output,
                    key_cols=sort_columns,
                    col_order=col_order,
                    chunksize=cfg.io.csv_chunksize,
                    sort_chunksize=cfg.io.csv_chunksize,
                    sep=cfg.io.csv_sep,
                    encoding=cfg.io.csv_encoding,
                    cfg=cfg,
                )
                logger.info("write_done", rows=rows_kept, path=str(csv_path))
            except OSError as exc:
                logger.error(
                    "write_fail",
                    error=str(exc),
                    path=str(output),
                )
                return 1

        if limit is not None:
            logger.info("process_limit", limit=processed_ids)

        if missing_required_columns:
            logger.warning(
                "validation_skipped",
                missing_columns=sorted(missing_required_columns),
            )
            validation_enabled = False
        elif optional_cols - present_columns:
            logger.warning(
                "optional_columns_missing",
                columns=sorted(optional_cols - present_columns),
            )

        if total_failures:
            errors.save(failure_path)

        assert csv_path is not None

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
            analyze_table_quality(csv_path, table_name=str(output.with_suffix("")))
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
        "ChEMBL activity data utilities",
        column="activity_id",
        chunk_size=5,
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
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
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
                "timeout": "activity.timeout",
                "column": "activity.column",
                "batch_size": "activity.batch_size",
                "limit": "activity.limit",
                "dry_run": "activity.dry_run",
            },
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
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
