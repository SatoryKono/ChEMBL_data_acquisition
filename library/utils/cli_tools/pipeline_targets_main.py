"""CLI wrapper for :mod:`library.pipeline_targets`.

The interface honours ``--base-path``, ``--input-dir`` and ``--output-dir`` when
resolving file locations.
"""

from __future__ import annotations

# ruff: noqa: E402
import argparse
import sys
from collections.abc import Iterable, Iterator, Sequence
from itertools import islice
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
    path_argument,
)
from library.config import Config, ConfigError, ensure_dirs, print_config
from library.io.paths import default_output_path
from library.io.readers import read_ids
from library.io.writers import write_csv
from library.log import logger
from library.pipeline_metadata import pipeline_metadata
from library.pipeline_targets import run_pipeline
from library.schemas import TargetsSchema
from library.schemas.targets import TARGETS_COLUMN_ORDER


class PipelineConfig(BaseModel):
    """Validated parameters for the target pipeline CLI."""

    model_config = ConfigDict(extra="forbid")

    input_csv: Path
    output_csv: Path | None = None
    final_out: Path | None = None
    raw_out: Path | None = None
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
        description="Run the target acquisition pipeline",
        parents=[root],
        conflict_handler="resolve",
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
    parser.add_argument(
        "--raw-out",
        dest="raw_out",
        type=path_argument,
        default=None,
        help="Optional CSV capturing the raw concatenated payload",
    )
    parser.add_argument(
        "--final-out",
        dest="final_out",
        type=path_argument,
        default=None,
        help="Destination CSV aligned to TargetsSchema",
    )
    parser.add_argument(
        "--out",
        dest="final_out",
        type=path_argument,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
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


def _frame_iterator(data: pd.DataFrame | Iterable[pd.DataFrame]) -> Iterator[pd.DataFrame]:
    """Return an iterator over ``data`` irrespective of its concrete type."""

    if isinstance(data, pd.DataFrame):
        yield data
        return
    yield from data


_PIPELINE_METADATA_COLUMNS = tuple(pipeline_metadata().keys())
_PLACEHOLDER = "-"


def _raw_dump(
    frames: Iterable[pd.DataFrame],
    path: Path,
    *,
    cfg: Config,
    sep: str | None = None,
    encoding: str | None = None,
) -> Path:
    """Persist ``frames`` to ``path`` preserving columns and row order."""

    sep = sep or cfg.io.csv_sep
    encoding = encoding or cfg.io.csv_encoding
    path.parent.mkdir(parents=True, exist_ok=True)
    iterator = iter(frames)
    try:
        first = next(iterator)
    except StopIteration:
        first = pd.DataFrame()
    columns = list(first.columns)
    first.to_csv(path, index=False, sep=sep, encoding=encoding)
    for chunk in iterator:
        if columns:
            chunk = chunk.reindex(columns=columns)
        else:
            columns = list(chunk.columns)
        chunk.to_csv(
            path,
            index=False,
            sep=sep,
            encoding=encoding,
            mode="a",
            header=False,
    )
    return path


class _RawStreamWriter:
    """Incrementally persist raw output chunks to disk."""

    def __init__(
        self,
        path: Path,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
    ) -> None:
        self.path = path
        self.sep = sep or cfg.io.csv_sep
        self.encoding = encoding or cfg.io.csv_encoding
        self._columns: list[str] | None = None
        self._written = False
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def write(self, chunk: pd.DataFrame) -> None:
        if self._columns is None:
            self._columns = list(chunk.columns)
            chunk.to_csv(
                self.path,
                index=False,
                sep=self.sep,
                encoding=self.encoding,
            )
            self._written = True
            return
        if self._columns:
            chunk = chunk.reindex(columns=self._columns)
        else:
            self._columns = list(chunk.columns)
        chunk.to_csv(
            self.path,
            index=False,
            sep=self.sep,
            encoding=self.encoding,
            mode="a",
            header=False,
        )
        self._written = True

    def close(self) -> Path:
        if not self._written:
            empty = pd.DataFrame(columns=self._columns or None)
            empty.to_csv(
                self.path,
                index=False,
                sep=self.sep,
                encoding=self.encoding,
            )
            self._written = True
        return self.path


def _replace_id_placeholders(
    df: pd.DataFrame, placeholder: str = _PLACEHOLDER
) -> pd.DataFrame:
    """Return ``df`` with placeholder values removed from identifier columns."""

    if df.empty:
        return df.copy()
    result = df.copy()
    id_columns = [column for column in result.columns if "id" in column.lower()]
    for column in id_columns:
        series = result[column]
        mask = series == placeholder
        if mask.any():
            result.loc[mask, column] = ""
    return result


def _prepare_final_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Align ``df`` to :data:`TARGETS_COLUMN_ORDER` with placeholder fills."""

    required_columns = [
        name for name, column in TargetsSchema.columns.items() if column.required
    ]
    if required_columns:
        required_schema = TargetsSchema.select_columns(required_columns)
        required_schema.validate(df, lazy=True)

    prepared = df.reindex(columns=TARGETS_COLUMN_ORDER, fill_value=_PLACEHOLDER)
    prepared = _replace_id_placeholders(prepared)
    return prepared


def _validate_and_write(
    frames: pd.DataFrame | Iterable[pd.DataFrame],
    path: Path,
    *,
    cfg: Config,
    sep: str | None = None,
    encoding: str | None = None,
) -> Path:
    """Validate ``df`` against :data:`TargetsSchema` and serialise to CSV."""

    def _validated_stream() -> Iterator[pd.DataFrame]:
        if not isinstance(frames, pd.DataFrame):
            iterator = iter(frames)
            yielded = False
            for chunk in iterator:
                yielded = True
                prepared = _prepare_final_frame(chunk)
                yield TargetsSchema.validate(prepared, lazy=True)
            if not yielded:
                prepared = _prepare_final_frame(pd.DataFrame())
                yield TargetsSchema.validate(prepared, lazy=True)
            return

        prepared = _prepare_final_frame(frames)
        yield TargetsSchema.validate(prepared, lazy=True)

    path.parent.mkdir(parents=True, exist_ok=True)
    write_csv(
        _validated_stream(),
        path,
        cfg=cfg,
        sep=sep,
        encoding=encoding,
        key_cols=["target_chembl_id"],
        col_order=TARGETS_COLUMN_ORDER,
    )
    return path


def _resolve_optional_output(
    path: Path | None | object,
    *,
    base_path: Path | None,
    output_dir: Path | None,
) -> Path | None:
    """Return an absolute path for optional outputs respecting CLI roots."""

    if path is None or path is argparse.SUPPRESS:
        return None
    candidate = path if isinstance(path, Path) else Path(str(path))
    if candidate.is_absolute():
        return candidate
    if output_dir is not None:
        return (output_dir / candidate).resolve()
    if base_path is not None:
        return (base_path / candidate).resolve()
    return candidate.resolve()


def _cached_chembl_fetch(
    chunks: Iterator[Iterable[str]],
    cfg: Config,
    *,
    chunk_size: int = 100,
    batch_size: int | None = None,
) -> Iterator[pd.DataFrame]:
    """Yield deterministic ChEMBL frames for the provided ``chunks``.

    The wrapper mimics the behaviour of the production fetcher which caches
    results on disk.  In this trimmed-down CLI variant we now expose a
    streaming interface so callers can process chunks incrementally without
    materialising the entire dataset at once.
    """

    def _stream() -> Iterator[pd.DataFrame]:
        yielded = False
        for frame in _chembl_frame_stream(chunks):
            yielded = True
            yield frame
        if not yielded:
            yield pd.DataFrame(
                {
                    "target_chembl_id": pd.Series(dtype="string"),
                    "source": pd.Series(dtype="string"),
                }
            )

    return _stream()


def run(cfg: Config, options: PipelineConfig) -> int:
    chunk_factory = lambda: _chunk_iterator(cfg, options)
    batch_size = options.batch_size if options.batch_size is not None else 100
    result = run_pipeline(
        chunk_factory,
        cfg,
        chembl_fetcher=_cached_chembl_fetch,
        chembl_kwargs={"chunk_size": options.chunk_size},
        batch_size=batch_size,
    )

    final_path = options.final_out or options.output_csv
    if final_path is None:
        final_path = default_output_path(options.input_csv, cfg.io)

    frame_iter = _frame_iterator(result.chembl)
    raw_writer: _RawStreamWriter | None = None
    if options.raw_out is not None:
        raw_writer = _RawStreamWriter(
            options.raw_out,
            cfg=cfg,
            sep=options.sep,
            encoding=options.encoding,
        )

    def _final_stream() -> Iterator[pd.DataFrame]:
        try:
            for chunk in frame_iter:
                if raw_writer is not None:
                    raw_chunk = chunk.drop(
                        columns=_PIPELINE_METADATA_COLUMNS, errors="ignore"
                    )
                    raw_writer.write(raw_chunk)
                yield chunk
        finally:
            if raw_writer is not None:
                raw_writer.close()

    _validate_and_write(
        _final_stream(),
        final_path,
        cfg=cfg,
        sep=options.sep,
        encoding=options.encoding,
    )
    logger.info("write_done", extra={"path": str(final_path)})
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    input_path = getattr(args, "input_csv", None)
    output_stem = Path(input_path).stem if input_path else None
    cli.prepare_io_paths(args, output_stem=output_stem)
    base_path = getattr(args, "base_path", None)
    output_dir = getattr(args, "output_dir", None)
    raw_candidate = getattr(args, "raw_out", None)
    resolved_raw = _resolve_optional_output(
        raw_candidate,
        base_path=base_path,
        output_dir=output_dir,
    )
    setattr(args, "raw_out", resolved_raw)
    final_candidate = getattr(args, "final_out", None)
    if final_candidate in (None, argparse.SUPPRESS):
        final_candidate = getattr(args, "output_csv", None)
    resolved_final = _resolve_optional_output(
        final_candidate,
        base_path=base_path,
        output_dir=output_dir,
    )
    setattr(args, "final_out", resolved_final)
    log_cfg.level = args.log_level
    log_cfg = LoggerConfig(level=log_cfg.level, run_id=log_cfg.run_id)
    logger_inst = configure_logger(log_cfg)
    pipeline_logger = logger_inst.bind(stage="pipeline")
    pipeline_logger.info("pipeline_start")
    try:
        cfg = cli.apply_config_overrides(args, parser, args.config)
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        pipeline_logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        pipeline_logger.info("pipeline_done", exit_code=1)
        return 1

    if args.print_config:
        print_config(cfg)
        return 0

    try:
        ensure_dirs(cfg)
        options = PipelineConfig.model_validate(
            {
                "input_csv": args.input_csv,
                "output_csv": args.output_csv,
                "final_out": args.final_out,
                "raw_out": args.raw_out,
                "column": args.column,
                "chunk_size": args.chunk_size,
                "batch_size": args.batch_size,
                "sep": args.sep,
                "encoding": args.encoding,
                "na_markers": args.na_markers,
                "limit": args.limit,
            }
        )
    except (ValueError, TypeError) as exc:
        pipeline_logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        pipeline_logger.info("pipeline_done", exit_code=1)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        pipeline_logger.error(
            "directory_setup_failed",
            error=str(exc),
        )
        pipeline_logger.info("pipeline_done", exit_code=1)
        return 1
    except ValidationError as exc:
        parser.error(str(exc))
    exit_code = run(cfg, options)
    pipeline_logger.info("pipeline_done", exit_code=exit_code)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
