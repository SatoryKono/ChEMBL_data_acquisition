"""Command line interface for retrieving ChEMBL activity data.

The module exposes a ``main`` entry point compatible with setuptools console
scripts as well as helpers that can be invoked directly from other
applications or tests.
"""

from __future__ import annotations

import argparse
import sys

from collections.abc import Iterable, Iterator, Sequence

from functools import partial
from itertools import islice
from pathlib import Path

try:
    from library.utils.bootstrap import ensure_project_root
except ModuleNotFoundError:  # pragma: no cover - fallback for direct execution
    def ensure_project_root() -> None:
        """Add the repository root to ``sys.path`` when executed as a script."""

        project_root = str(Path(__file__).resolve().parent.parent)
        if project_root not in sys.path:
            sys.path.insert(0, project_root)

if __package__ in {None, ""}:
    ensure_project_root()

import pandas as pd
import requests

from library import chembl_library as cl
from library import cli
from library import io
from library.clients import ChemblClient
from library.csv_utils import write_csv_chunks_deterministic  # re-exported for tests
from library.pipeline_helpers import (
    ChunkedFetchConfig,
    CsvWriterConfig,
    prepare_chunked_pipeline,
)
from library.rate_limiter import get_global_limiter

from library.cli import (
    LoggerConfig,
    positive_int,
)
from library.cli import (
    build_parser as base_parser,
)
from library.cli_utils import (
    PipelineError,
    resolve_invocation,
    run_cli_command,
    run_pipeline,
)
from library.cli_utils import (
    file_sha256 as _cli_file_sha256,
)
from library.cli_utils import (
    write_meta_yaml as _cli_write_meta_yaml,
)
from library.config import Config, _serialize_paths
from library.log import logger
from library.pipeline_metadata import add_pipeline_metadata
from library.processing.activity import (
    apply_activity_annotations,
    compute_activity_bounds,
)
from library.rate_limiter import get_global_limiter
from library.table_quality import analyze_table_quality
from library.validation import validate_activities
from schemas import ActivitiesSchema, configure_activity_schema, normalize_activities

DEFAULT_INPUT_NAME = "activity.csv"
DEFAULT_OUTPUT_STEM = "activities"
PROGRAM_NAME = Path(__file__).with_suffix("").name


def _args_invocation(args: argparse.Namespace) -> tuple[str, ...]:
    invocation = getattr(args, "invocation", None)
    if invocation is None:
        return (PROGRAM_NAME,)
    return tuple(str(arg) for arg in invocation)


file_sha256 = _cli_file_sha256
write_meta_yaml = _cli_write_meta_yaml
configure_logger = cli.configure_logger

__all__ = (
    "file_sha256",
    "write_meta_yaml",
    "configure_logger",
)


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute activity retrieval from the ChEMBL API.

    The resulting CSV places columns defined in :data:`ActivitiesSchema`
    first, preserving their declared order. Any additional fields appear
    afterwards sorted alphabetically.

    Parameters
    ----------
    cfg : Config
        Application configuration providing API credentials, retry strategy
        and CSV export options.
    args : argparse.Namespace
        Parsed command-line arguments. ``args.limit`` constrains the number of
        identifiers processed, while ``args.dry_run`` skips network calls and
        file generation.

    Returns
    -------
    int
        ``0`` on success, non-zero when validation or I/O failures are
        encountered. Upstream API errors are logged and converted into a
        failure code by :func:`library.cli_utils.run_pipeline`.
    """
    limit = cfg.activity.limit
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="activity.limit", limit=limit)
        return 1

    if cfg.activity.dry_run:
        expected = limit if limit is not None else 0
        logger.info("dry_run", limit=expected)
        return 0

    try:
        ids_iter = io.read_ids(args.input_csv, column=cfg.activity.column, cfg=cfg.io)
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

    processed_ids = 0

    def _iter_ids() -> Iterator[str]:
        nonlocal processed_ids
        for identifier in ids_iter:
            processed_ids += 1
            yield identifier

    limited_ids: Iterator[str]
    if limit is not None:
        limited_ids = islice(_iter_ids(), limit)
    else:
        limited_ids = _iter_ids()

    enrichment_cfg = cfg.activity_enrichment
    extra_columns: list[str] = []
    action_cfg = enrichment_cfg.action_type
    configure_activity_schema(action_cfg.metrics)
    if action_cfg.enabled or action_cfg.log_missing or action_cfg.log_distribution:
        extra_columns.append(action_cfg.column)
    extra_kwargs = {"extra_columns": extra_columns} if extra_columns else {}

    invocation = _args_invocation(args)

    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    failure_path = Path(output).with_name(f"{Path(output).stem}_failure_cases.csv")


    def _compute_bounds(frame: pd.DataFrame) -> pd.DataFrame:
        return compute_activity_bounds(frame, cfg.activity_bounds)

    def _apply_annotations(frame: pd.DataFrame) -> pd.DataFrame:
        return apply_activity_annotations(
            frame,
            action_cfg=enrichment_cfg.action_type,
            properties_cfg=enrichment_cfg.activity_properties,
        )

    metadata_hooks = [
        normalize_activities,
        add_pipeline_metadata,
        _compute_bounds,
        _apply_annotations,
    ]

    validators = [partial(validate_activities, return_result=True)]



    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: Sequence[str],
        key_cols: Sequence[str],
    ) -> Path:
        sort_columns = list(key_cols) or sorted(col_order)
        output_path = write_csv_chunks_deterministic(
            chunks,
            destination,
            key_cols=sort_columns,
            col_order=list(col_order),
            chunksize=cfg.io.csv_chunksize,
            sort_chunksize=cfg.io.csv_chunksize,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
        )
        path_obj = Path(output_path)
        if not path_obj.exists():
            raise PipelineError(f"writer failed to create output: {path_obj}")
        return path_obj



    table_quality = partial(
        analyze_table_quality,
        table_name=str(Path(output).with_suffix("")),
        destination_dir=Path(output).parent,
    )

    global_limiter = get_global_limiter(
        cfg.rate.global_rps,
        cfg.rate.global_burst,
    )

    with ChemblClient(
        cfg.api,
        cfg.retry,
        cfg.chembl,
        global_limiter=global_limiter,
    ) as client:

        def fetch_chunk(chunk_ids: Sequence[str]) -> pd.DataFrame:
            try:
                return cl.get_activities(
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
                    extra={
                        "msg": str(exc),
                        "chunk_ids": list(chunk_ids),
                    },
                    error=str(exc),
                    batch_size=cfg.activity.batch_size,
                    timeout=cfg.activity.timeout,
                )
                raise PipelineError(str(exc)) from exc

        worker_count = getattr(cfg.activity, "workers", 1) or 1
        fetch_config = ChunkedFetchConfig(
            ids=limited_ids,
            chunk_size=cfg.activity.batch_size,
            workers=max(1, worker_count),
        )

        def _write_chunks(
            chunks: Iterable[pd.DataFrame],
            destination: Path,
            *,
            col_order: Sequence[str] | None,
            key_cols: Sequence[str] | None,
            **kwargs,
        ) -> Path:
            path = write_csv_chunks_deterministic(
                chunks,
                destination,
                col_order=col_order,
                key_cols=key_cols,
                **kwargs,
            )
            path_obj = Path(path)
            if not path_obj.exists():
                raise PipelineError(f"writer failed to create output: {path_obj}")
            return path_obj

        writer_config = CsvWriterConfig(
            writer=writer,
            kwargs={},
            ensure_destination=False,

        )

        fetcher, writer = prepare_chunked_pipeline(
            fetch_config=fetch_config,
            fetch_chunk=fetch_chunk,
            csv_writer=writer_config,
        )

        exit_code = run_pipeline(
            fetcher=fetcher,
            schema=ActivitiesSchema,
            schema_name="ActivitiesSchema",
            validators=validators,
            metadata_hooks=metadata_hooks,
            writer=writer,
            output_path=output,
            failure_path=failure_path,
            command=" ".join(_args_invocation(args)),

            config_snapshot=_serialize_paths(cfg.to_dict()),
            inputs={"input_csv": str(args.input_csv)},
            key_columns=["activity_id"],
            table_quality=table_quality,
            cfg=cfg,
            logger=logger,
        )

    if limit is not None:
        logger.info("process_limit", limit=processed_ids)

    return exit_code


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the activity pipeline handling ``--skip-existing`` semantics."""

    output_path = Path(
        args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    )
    args.output_csv = output_path
    if args.skip_existing and output_path.exists() and not args.force:
        logger.info("pipeline_skip_existing", output=str(output_path))
        return 0
    return run_chembl(cfg, args)


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        A tuple containing the fully configured parser and the logging
        configuration populated with defaults.
    """
    parser, log_cfg = base_parser(
        "ChEMBL activity data utilities",
        column="activity_id",
        chunk_size=5,
        size_option="--batch-size",
        size_dest="batch_size",
    )
    parser.prog = PROGRAM_NAME
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
        help=(
            "Maximum number of identifiers to process; use 0 to skip processing"
        ),
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
    parser.add_argument(
        "--workers",
        type=positive_int,
        default=1,
        help="Number of worker threads fetching activities in parallel",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the activity pipeline with optional argument overrides.

    Parameters
    ----------
    argv : Sequence[str] | None, optional
        Command-line arguments to parse. When ``None`` the values from
        :data:`sys.argv` are used.

    Returns
    -------
    int
        ``0`` on success, non-zero otherwise. Errors are logged with
        structured context for easier diagnosis.

    Raises
    ------
    SystemExit
        Raised when argument parsing fails, mirroring ``argparse`` behaviour.
    """
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    args.invocation = resolve_invocation(parser.prog, argv)
    cli.prepare_io_paths(
        args,
        input_default=DEFAULT_INPUT_NAME,
        output_stem=DEFAULT_OUTPUT_STEM,
    )
    if args.limit == 0:
        logger.info("pipeline_skip_limit", limit=args.limit)
        return 0
    if args.limit is not None and args.limit < 0:
        # Reject negative limits early to provide clear CLI feedback.
        parser.error("--limit must be zero or a positive integer")
    if args.offset < 0:
        parser.error("--offset must be zero or a positive integer")
    return run_cli_command(
        args=args,
        parser=parser,
        log_cfg=log_cfg,
        mapping={
            "timeout": "activity.timeout",
            "column": "activity.column",
            "batch_size": "activity.batch_size",
            "limit": "activity.limit",
            "offset": "activity.offset",
            "dry_run": "activity.dry_run",
            "workers": "activity.workers",
        },
        run=run,
        logger=logger,
    )


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
