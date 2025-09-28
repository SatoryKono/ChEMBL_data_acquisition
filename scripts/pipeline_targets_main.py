"""CLI wrapper for :mod:`library.pipeline_targets`."""

from __future__ import annotations

# ruff: noqa: E402
import argparse
import sys
from collections.abc import Iterable, Iterator, Sequence
from itertools import islice
from pathlib import Path

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, ValidationError

from library import cli
from library.chembl_client import _chunked
from library.cli import (
    LoggerConfig,
    build_root_parser,
    configure_logger,
)
from library.config import Config, ensure_dirs, print_config
from library.io import default_output_path, read_ids, write_csv
from library.log import logger
from library.pipeline_metadata import add_pipeline_metadata
from library.pipeline_targets import PipelineResult, run_pipeline


class PipelineConfig(BaseModel):
    """Validated parameters for the target pipeline CLI."""

    model_config = ConfigDict(extra="forbid")

    input_csv: Path
    output_csv: Path | None = None
    column: str = Field("target_chembl_id")
    chunk_size: int = Field(100, ge=1)
    batch_size: int | None = Field(100, ge=1)
    sep: str | None = None
    encoding: str | None = None
    na_markers: list[str] | None = None
    limit: int | None = Field(default=None, ge=1)


def _positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:  # pragma: no cover - delegated to argparse
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def _optional_batch(value: str) -> int | None:
    if value.lower() in {"none", "null"}:
        return None
    return _positive_int(value)


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the argument parser for the pipeline CLI."""

    root, _shared, log_cfg = build_root_parser()
    parser = argparse.ArgumentParser(
        description="Run the target acquisition pipeline", parents=[root]
    )
    parser.add_argument(
        "--column",
        default="target_chembl_id",
        help="Column in the input CSV containing target identifiers",
    )
    parser.add_argument(
        "--chunk-size",
        type=_positive_int,
        default=100,
        help="Number of identifiers grouped per chunk",
    )
    parser.add_argument(
        "--batch-size",
        type=_optional_batch,
        default=100,
        help="Batch size forwarded to fetchers supporting it",
    )
    parser.add_argument(
        "--limit",
        type=_positive_int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    parser.add_argument(
        "--na-marker",
        action="append",
        dest="na_markers",
        help="Additional value interpreted as missing",
    )
    return parser, log_cfg


def _chunk_iterator(cfg: Config, options: PipelineConfig) -> Iterator[Iterable[str]]:
    ids = read_ids(
        options.input_csv,
        column=options.column,
        cfg=cfg.io,
        sep=options.sep,
        encoding=options.encoding,
        na_markers=options.na_markers,
    )
    if options.limit is not None:
        ids = islice(ids, options.limit)
    yield from _chunked(ids, options.chunk_size)


def _cached_chembl_fetch(
    chunks: Iterator[Iterable[str]],
    cfg: Config,
    *,
    chunk_size: int = 100,
    batch_size: int | None = None,
) -> pd.DataFrame:
    """Return a deterministic frame for the provided ``chunks``.

    The wrapper mimics the behaviour of the production fetcher which caches
    results on disk.  In this trimmed-down CLI variant we simply flatten the
    incoming chunk iterator into a dataframe while keeping the API compatible
    with :func:`library.pipeline_targets.run_pipeline`.
    """

    ids = [item for chunk in chunks for item in chunk]
    df = pd.DataFrame({"target_chembl_id": ids})
    if df.empty:
        return df
    df["source"] = "chembl"
    return df


def _write_outputs(
    cfg: Config, options: PipelineConfig, result: PipelineResult
) -> Path:
    output = options.output_csv or default_output_path(options.input_csv, cfg.io)
    annotated = add_pipeline_metadata(result.chembl)
    write_csv(
        annotated,
        output,
        cfg=cfg,
        sep=options.sep,
        encoding=options.encoding,
    )
    return output


def run(cfg: Config, options: PipelineConfig) -> int:
    result = run_pipeline(
        lambda: _chunk_iterator(cfg, options),
        cfg,
        chembl_fetcher=_cached_chembl_fetch,
        chembl_kwargs={"chunk_size": options.chunk_size},
        batch_size=options.batch_size,
    )
    output = _write_outputs(cfg, options, result)
    logger.info("write_done", extra={"path": str(output)})
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    log_cfg.level = args.log_level
    log_cfg = LoggerConfig(level=log_cfg.level, run_id=log_cfg.run_id)
    logger_inst = configure_logger(log_cfg)
    pipeline_logger = logger_inst.bind(stage="pipeline")
    pipeline_logger.info("pipeline_start")
    cfg = cli.apply_config_overrides(args, parser, args.config)
    ensure_dirs(cfg)
    print_config(cfg)
    if args.print_config:
        return 0
    try:
        options = PipelineConfig.model_validate(
            {
                "input_csv": args.input_csv,
                "output_csv": args.output_csv,
                "column": args.column,
                "chunk_size": args.chunk_size,
                "batch_size": args.batch_size,
                "sep": args.sep,
                "encoding": args.encoding,
                "na_markers": args.na_markers,
                "limit": args.limit,
            }
        )
    except ValidationError as exc:
        parser.error(str(exc))
    exit_code = run(cfg, options)
    pipeline_logger.info("pipeline_done", exit_code=exit_code)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
