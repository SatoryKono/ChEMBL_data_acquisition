"""CLI wrapper for :mod:`library.pipeline_targets`.

The interface honours ``--base-path``, ``--input-dir`` and ``--output-dir`` when
resolving file locations.
"""

from __future__ import annotations

# ruff: noqa: E402
import argparse
import sys
from collections.abc import Iterable, Iterator, Sequence
from itertools import chain, islice
from pathlib import Path

if __package__ in {None, ""}:
    project_root = Path(__file__).resolve().parents[3]
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

from library.utils import bootstrap

bootstrap.ensure_project_root()

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, ValidationError

from library import cli
from library.clients import _chunked
from library.cli import (
    LoggerConfig,
    build_root_parser,
    configure_logger,
)
from library.config import Config, ensure_dirs, print_config
from library.io.paths import default_output_path
from library.io.readers import read_ids
from library.io.writers import write_csv
from library.log import logger
from library.pipeline_metadata import add_pipeline_metadata
from library.pipeline_targets import run_pipeline


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
    limit: int | None = Field(default=None, ge=0)


def _non_negative_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:  # pragma: no cover - delegated to argparse
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return parsed


def _positive_int(value: str) -> int:
    parsed = _non_negative_int(value)
    if parsed == 0:
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
        type=_non_negative_int,
        default=None,
        help="Maximum number of identifiers to process (0 skips processing)",
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


def _chembl_frame_stream(chunks: Iterator[Iterable[str]]) -> Iterator[pd.DataFrame]:
    """Yield one dataframe per ``chunks`` entry with ChEMBL metadata added."""

    for chunk in chunks:
        ids = list(chunk)
        if not ids:
            continue
        frame = pd.DataFrame({"target_chembl_id": ids})
        frame["source"] = "chembl"
        yield frame


def _cached_chembl_fetch(
    chunks: Iterator[Iterable[str]],
    cfg: Config,
    *,
    chunk_size: int = 100,
    batch_size: int | None = None,
) -> pd.DataFrame:
    """Return a deterministic frame for the provided ``chunks``.

    The wrapper mimics the behaviour of the production fetcher which caches
    results on disk.  In this trimmed-down CLI variant we lazily convert the
    incoming chunk iterator into dataframes and concatenate them on demand
    while keeping the API compatible with
    :func:`library.pipeline_targets.run_pipeline`.
    """

    frames = _chembl_frame_stream(chunks)
    try:
        first = next(frames)
    except StopIteration:
        return pd.DataFrame(
            {
                "target_chembl_id": pd.Series(dtype="string"),
                "source": pd.Series(dtype="string"),
            }
        )
    return pd.concat(chain([first], frames), ignore_index=True)


def _write_outputs(
    cfg: Config, options: PipelineConfig, frames: Iterable[pd.DataFrame]
) -> Path:
    output = options.output_csv or default_output_path(options.input_csv, cfg.io)
    iterator = iter(frames)
    try:
        first = next(iterator)
    except StopIteration:
        empty = add_pipeline_metadata(
            pd.DataFrame(
                {
                    "target_chembl_id": pd.Series(dtype="string"),
                    "source": pd.Series(dtype="string"),
                }
            )
        )
        write_csv(
            empty,
            output,
            cfg=cfg,
            sep=options.sep,
            encoding=options.encoding,
        )
        return output

    annotated_first = add_pipeline_metadata(first)
    annotated_rest = (add_pipeline_metadata(chunk) for chunk in iterator)
    annotated_chunks = chain([annotated_first], annotated_rest)
    write_csv(
        annotated_chunks,
        output,
        cfg=cfg,
        sep=options.sep,
        encoding=options.encoding,
    )
    return output


def run(cfg: Config, options: PipelineConfig) -> int:
    chunk_factory = lambda: _chunk_iterator(cfg, options)
    _ = run_pipeline(
        chunk_factory,
        cfg,
        chembl_fetcher=_cached_chembl_fetch,
        chembl_kwargs={"chunk_size": options.chunk_size},
        batch_size=options.batch_size,
    )
    output = _write_outputs(cfg, options, _chembl_frame_stream(chunk_factory()))
    logger.info("write_done", extra={"path": str(output)})
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    input_path = getattr(args, "input_csv", None)
    output_stem = Path(input_path).stem if input_path else None
    cli.prepare_io_paths(args, output_stem=output_stem)
    log_cfg.level = args.log_level
    log_cfg = LoggerConfig(level=log_cfg.level, run_id=log_cfg.run_id)
    logger_inst = configure_logger(log_cfg)
    pipeline_logger = logger_inst.bind(stage="pipeline")
    pipeline_logger.info("pipeline_start")
    cfg = cli.apply_config_overrides(args, parser, args.config)
    if args.print_config:
        print_config(cfg)
        return 0
    ensure_dirs(cfg)
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
