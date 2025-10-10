"""CLI wrapper for :mod:`library.pipelines.target.pipeline`.

The interface honours ``--base-path``, ``--input-dir`` and ``--output-dir`` when
resolving file locations.
"""

from __future__ import annotations

# ruff: noqa: E402
import argparse
import sys
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import replace
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
from library.cli import (
    LoggerConfig,
    build_root_parser,
    configure_logger,
    path_argument,
)
from library.cli.run_context import compute_generated_at
from library.clients import _chunked
from library.common.log import logger
from library.config import Config, ConfigError, ensure_dirs, print_config
from library.io import write_meta_yaml
from library.io.paths import default_output_path
from library.io.readers import read_ids
from library.io.writers import write_csv
from library.pipelines.common import pipeline_metadata
from library.pipelines.target.pipeline import run_pipeline
from library.postprocessing import target as target_pp
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
    raw_format: str = Field("csv")
    id_cols: list[str] | None = None
    no_reindex_raw: bool = False
    normalize_at_export: bool = True


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
        "--raw-format",
        dest="raw_format",
        choices=("csv", "parquet"),
        default="csv",
        help="Format used when writing raw payload exports",
    )
    parser.add_argument(
        "--id-cols",
        dest="id_cols",
        nargs="+",
        default=None,
        help="Identifier columns used for deterministic ordering",
    )
    parser.add_argument(
        "--no-reindex-raw",
        dest="no_reindex_raw",
        action="store_true",
        default=False,
        help="Skip column reindexing when exporting the raw dataset",
    )
    parser.add_argument(
        "--normalize-at-export",
        "--no-normalize-at-export",
        dest="normalize_at_export",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Apply schema alignment before writing the final output. "
            "Use --no-normalize-at-export to persist the raw payload."
        ),
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


def _frame_iterator(
    data: pd.DataFrame | Iterable[pd.DataFrame],
) -> Iterator[pd.DataFrame]:
    """Return an iterator over ``data`` irrespective of its concrete type."""

    if isinstance(data, pd.DataFrame):
        yield data
        return
    yield from data


_PIPELINE_METADATA_COLUMNS = tuple(pipeline_metadata().keys())
_PLACEHOLDER = "-"


def _raw_dump(
    frames: Iterable[pd.DataFrame] | pd.DataFrame,
    path: Path,
    *,
    cfg: Config,
    sep: str | None = None,
    encoding: str | None = None,
    raw_format: str = "csv",
    reindex_columns: bool = True,
) -> Path:
    """Persist ``frames`` to ``path`` preserving columns and row order."""

    writer = PipelineTargetsRawWriter(
        path,
        cfg=cfg,
        sep=sep,
        encoding=encoding,
        raw_format=raw_format,
        reindex_columns=reindex_columns,
    )
    if isinstance(frames, pd.DataFrame):
        writer.write(frames)
    else:
        for chunk in frames:
            writer.write(chunk)
    return writer.close()


class PipelineTargetsRawWriter:
    """Incrementally persist raw output chunks to disk."""

    def __init__(
        self,
        path: Path,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
        raw_format: str = "csv",
        reindex_columns: bool = True,
    ) -> None:
        self.path = path
        self._cfg = cfg
        self.sep = sep or cfg.io.csv_sep
        self.encoding = encoding or cfg.io.csv_encoding
        normalized_format = (raw_format or "csv").lower()
        if normalized_format not in {"csv", "parquet"}:
            logger.warning(
                "unsupported_raw_format",
                format=raw_format,
                fallback="csv",
            )
            normalized_format = "csv"
        self.raw_format = normalized_format
        self._reindex = reindex_columns
        self._columns: list[str] | None = None
        self._dtypes: dict[str, str] | None = None
        self._rows_written = 0
        self._frames: list[pd.DataFrame] | None = (
            [] if self.raw_format == "parquet" else None
        )
        self._destination_opened = False
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def write(self, chunk: pd.DataFrame) -> None:
        if chunk.empty and self._rows_written == 0 and self._columns is None:
            if self._reindex:
                self._columns = []
            else:
                self._columns = list(chunk.columns)
            return

        if self._columns is None:
            self._columns = list(chunk.columns)
        else:
            new_columns = [c for c in chunk.columns if c not in self._columns]
            if new_columns:
                if self.raw_format == "parquet" or self._rows_written == 0:
                    merged = list(self._columns)
                    merged.extend(col for col in new_columns if col not in merged)
                    self._columns = merged
                else:
                    raise OSError("raw_dump_inconsistent_columns")

        working = chunk
        if self._columns is not None:
            working = working.reindex(columns=self._columns)

        if self.raw_format == "parquet":
            if self._frames is None:
                self._frames = []
            self._frames.append(working.copy())
        else:
            mode = "w" if self._rows_written == 0 else "a"
            header = self._rows_written == 0
            working.to_csv(
                self.path,
                index=False,
                sep=self.sep,
                encoding=self.encoding,
                mode=mode,
                header=header,
                lineterminator="\n",
            )
            self._destination_opened = True

        if self._columns is not None:
            self._dtypes = {col: str(working[col].dtype) for col in self._columns}
        else:
            self._dtypes = {col: str(dtype) for col, dtype in working.dtypes.items()}

        self._rows_written += len(working)

    def close(self) -> Path:
        if self.raw_format == "parquet":
            frames = self._frames or []
            if frames:
                combined = pd.concat(frames, ignore_index=True)
            else:
                combined = pd.DataFrame(columns=self._columns or [])
            try:
                combined.to_parquet(self.path, index=False)
            except (ImportError, ValueError) as exc:
                raise OSError(f"failed_to_write_parquet: {exc}") from exc
        else:
            if not self._destination_opened:
                empty = pd.DataFrame(columns=self._columns or None)
                empty.to_csv(
                    self.path,
                    index=False,
                    sep=self.sep,
                    encoding=self.encoding,
                    lineterminator="\n",
                )
            columns = list(self._columns or [])
            dtypes = self._dtypes or {}
            write_meta_yaml(
                self.path,
                self._cfg,
                columns=columns,
                dtypes=dtypes or None,
            )
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
    id_cols: Sequence[str] | None = None,
    normalize: bool = True,
) -> Path:
    """Validate ``df`` against :data:`TargetsSchema` and serialise to CSV."""

    resolved_key_cols = list(id_cols) if id_cols else ["target_chembl_id"]

    def _ensure_key_columns(frame: pd.DataFrame) -> pd.DataFrame:
        if not resolved_key_cols:
            return frame
        missing = [col for col in resolved_key_cols if col not in frame.columns]
        if not missing:
            return frame
        ordered = list(frame.columns) + missing
        return frame.reindex(columns=ordered)

    def _validated_stream() -> Iterator[pd.DataFrame]:
        if not isinstance(frames, pd.DataFrame):
            iterator = iter(frames)
            yielded = False
            for chunk in iterator:
                yielded = True
                if normalize:
                    prepared = _prepare_final_frame(chunk)
                    yield TargetsSchema.validate(prepared, lazy=True)
                else:
                    yield _ensure_key_columns(chunk)
            if not yielded:
                if normalize:
                    prepared = _prepare_final_frame(pd.DataFrame())
                    yield TargetsSchema.validate(prepared, lazy=True)
                else:
                    yield pd.DataFrame(columns=resolved_key_cols)
            return

        if normalize:
            prepared = _prepare_final_frame(frames)
            yield TargetsSchema.validate(prepared, lazy=True)
        else:
            yield _ensure_key_columns(frames)

    path.parent.mkdir(parents=True, exist_ok=True)
    key_cols = resolved_key_cols
    col_order: Sequence[str] | None = TARGETS_COLUMN_ORDER if normalize else None
    write_csv(
        _validated_stream(),
        path,
        cfg=cfg,
        sep=sep,
        encoding=encoding,
        key_cols=key_cols,
        col_order=col_order,
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
    def chunk_factory() -> Iterator[pd.DataFrame]:
        return _chunk_iterator(cfg, options)

    batch_size = options.batch_size if options.batch_size is not None else 100
    raw_format = (options.raw_format or "csv").lower()
    if raw_format not in {"csv", "parquet"}:
        logger.warning("unsupported_raw_format", format=raw_format, fallback="csv")
        raw_format = "csv"
    reindex_raw = not options.no_reindex_raw
    id_cols = list(options.id_cols) if options.id_cols is not None else None
    normalize = options.normalize_at_export
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
    raw_writer: PipelineTargetsRawWriter | None = None
    if options.raw_out is not None:
        raw_writer = PipelineTargetsRawWriter(
            options.raw_out,
            cfg=cfg,
            sep=options.sep,
            encoding=options.encoding,
            raw_format=raw_format,
            reindex_columns=reindex_raw,
        )

    def _final_stream() -> Iterator[pd.DataFrame]:
        try:
            metadata_columns = list(_PIPELINE_METADATA_COLUMNS)
            for chunk in frame_iter:
                raw_chunk = chunk.drop(columns=metadata_columns, errors="ignore")
                if raw_writer is not None:
                    raw_writer.write(raw_chunk)
                yield chunk if normalize else raw_chunk
        finally:
            if raw_writer is not None:
                raw_writer.close()

    stream = _final_stream()
    skip_final_write = (
        not normalize
        and options.raw_out is not None
        and final_path == options.raw_out
        and raw_format == "parquet"
    )
    if skip_final_write:
        for _ in stream:
            pass
    else:
        _validate_and_write(
            stream,
            final_path,
            cfg=cfg,
            sep=options.sep,
            encoding=options.encoding,
            id_cols=id_cols,
            normalize=normalize,
        )
    if final_path.suffix.lower() == ".csv":
        target_pp.process_targets(str(final_path), verbose=True)
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
    args.raw_out = resolved_raw
    final_candidate = getattr(args, "final_out", None)
    if final_candidate in (None, argparse.SUPPRESS):
        final_candidate = getattr(args, "output_csv", None)
    resolved_final = _resolve_optional_output(
        final_candidate,
        base_path=base_path,
        output_dir=output_dir,
    )
    args.final_out = resolved_final
    run_id_value = getattr(args, "run_id", None)
    if isinstance(run_id_value, str):
        run_id_value = run_id_value.strip() or None
    if run_id_value is not None:
        updated_run_id = run_id_value
    else:
        updated_run_id = log_cfg.run_id
    resolved_run_id = str(updated_run_id)
    log_level_value = str(args.log_level)
    new_generated_at = compute_generated_at(
        date_token=None,
        run_id=resolved_run_id,
        seed_parts=("create_logger_config", log_level_value.upper()),
    )
    log_cfg = replace(
        log_cfg,
        level=log_level_value,
        run_id=resolved_run_id,
        generated_at=new_generated_at,
    )
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
                "raw_format": args.raw_format,
                "id_cols": args.id_cols,
                "no_reindex_raw": args.no_reindex_raw,
                "normalize_at_export": args.normalize_at_export,
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
