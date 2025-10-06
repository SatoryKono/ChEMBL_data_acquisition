"""Document pipeline execution helpers used by the CLI and tests."""

from __future__ import annotations

import argparse
import inspect
import sys
import tempfile
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from itertools import chain, islice
from numbers import Integral, Real
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from pandera.errors import SchemaErrors

from library import io
from library.clients import ChemblClient
from library.cli import ConfigMetadata
from library.cli.metadata import prepare_option
from library.common.csv_utils import write_csv_chunks_deterministic
from library.common.log import logger
from library.common.rate_limiter import get_global_limiter
from library.common.sidecar import SidecarErrors
from library.config import Config, _serialize_paths
from library.integration import chembl_library as cl
from library.pipelines.common import add_pipeline_metadata
from library.pipelines.document import postprocessing as dp
from library.pipelines.document.pipeline import (
    DOCUMENT_SCHEMA_COLUMNS,
    DocumentQualityAccumulator,
    build_dataframe,
    build_quality_report,
    dataframe_to_strings,
    merge_with_chembl,
    normalise_doi,
)
from library.pipelines.document.service import (
    DocumentPipeline,
    FallbackDoiMetrics,
    FallbackDoiState,
)
from library.postprocessing import document as document_export_postprocessing
from library.postprocessing.document import preprocess_documents_csv
from library.qa.table_quality import TableQualityProfiler
from library.reporting.run_manifest import (
    QualityAnalysisError,
    QualityReportError,
    finalise_csv_output,
)
from library.schemas import DocumentsSchema, normalize_documents
from library.schemas.document_spec import DOCUMENT_EXPORT_COLUMNS

_EXPORT_COLUMNS: tuple[str, ...] = DOCUMENT_EXPORT_COLUMNS


def _resolve_numeric_export_columns() -> tuple[str, ...]:
    """Return schema columns that should retain their numeric dtype."""

    numeric_columns: list[str] = []
    for name, column in DocumentsSchema.columns.items():
        dtype = getattr(column, "dtype", None)
        if dtype is None:
            continue
        dtype_name = str(dtype).lower()
        if dtype_name.startswith(("int", "float")):
            numeric_columns.append(name)
    return tuple(numeric_columns)


_NUMERIC_EXPORT_COLUMNS: tuple[str, ...] = _resolve_numeric_export_columns()

_EXPORT_COLUMN_RENAMES = {
    "document_chembl_id": "ChEMBL.document_chembl_id",
    "title": "ChEMBL.title",
    "abstract": "ChEMBL.abstract",
    "doi": "ChEMBL.doi",
    "year": "ChEMBL.year",
    "journal": "ChEMBL.journal",
    "journal_abbrev": "ChEMBL.journal_abbrev",
    "volume": "ChEMBL.volume",
    "issue": "ChEMBL.issue",
    "first_page": "ChEMBL.first_page",
    "last_page": "ChEMBL.last_page",
    "pubmed_id": "ChEMBL.pubmed_id",
    "authors": "ChEMBL.authors",
    "source": "ChEMBL.source",
    "publication_type_score_review": "publication_review_score",
    "publication_type_score_experimental": "publication_experimental_score",
}

_EXPORT_COALESCE_SOURCES = {
    "OpenAlex.PMID": ["OpenAlex.PMID", "PubMed.PMID", "scholar.PMID"],
    "OpenAlex.DOI": ["OpenAlex.DOI", "PubMed.DOI", "scholar.DOI", "doi_normalised"],
    "crossref.DOI": ["crossref.DOI", "doi_normalised", "PubMed.DOI", "scholar.DOI"],
}

_EXPORT_SORT_FALLBACK = [
    "ChEMBL.document_chembl_id",
    "PubMed.PMID",
    "scholar.PMID",
    "OpenAlex.PMID",
    "ChEMBL.pubmed_id",
]

_EXPORT_STREAM_CHUNK_SIZE = 10_000


def _coalesce_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.Series:
    """Return the first non-empty value across ``columns`` for each row."""

    result = pd.Series("", index=df.index, dtype=object)
    for col in columns:
        if col not in df.columns:
            continue
        values = df[col].fillna("").astype(str)
        mask = result.eq("")
        if mask.any():
            result.loc[mask] = values.loc[mask]
    return result


def _resolve_duplicate_column(frame: pd.DataFrame, column: str) -> pd.Series:
    """Return a single series for ``column`` consolidating duplicate columns."""

    if column not in frame.columns:
        return pd.Series(index=frame.index, dtype=object)

    selected = frame[column]
    if isinstance(selected, pd.Series):
        return selected

    if selected.shape[1] == 1:
        return selected.iloc[:, 0]

    consolidated = (
        selected.replace("", pd.NA)
        .bfill(axis=1)
        .ffill(axis=1)
        .iloc[:, 0]
        .fillna("")
    )
    return consolidated


def _collapse_duplicate_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Return ``frame`` with duplicate-named columns merged into single columns."""

    if frame.columns.empty:
        return frame.copy()

    column_order = list(dict.fromkeys(frame.columns))
    collapsed_columns: dict[Any, pd.Series] = {}
    for column in column_order:
        collapsed_columns[column] = _resolve_duplicate_column(frame, column)

    if not collapsed_columns:
        return frame.iloc[:, :0].copy()

    return pd.DataFrame(collapsed_columns, index=frame.index)


def _merge_preferred_series(target_series: pd.Series, source_series: pd.Series) -> pd.Series:
    """Return ``target_series`` with missing entries populated from ``source_series``."""

    if target_series.empty and source_series.empty:
        return target_series.copy()

    if not target_series.index.equals(source_series.index):
        source_series = source_series.reindex(target_series.index)

    combined = target_series.copy()

    missing_mask = target_series.isna()
    if missing_mask.any():
        combined.loc[missing_mask] = source_series.loc[missing_mask]

    if pd.api.types.is_object_dtype(target_series.dtype) or pd.api.types.is_string_dtype(
        target_series.dtype
    ):
        empty_mask = target_series.fillna("").eq("")
        if empty_mask.any():
            combined.loc[empty_mask] = source_series.loc[empty_mask]

    return combined


def _prepare_export_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Rename and project columns to match the export schema."""

    frame = _collapse_duplicate_columns(df.copy())

    rename_map: dict[str, str] = {}
    for source, target in _EXPORT_COLUMN_RENAMES.items():
        if source not in frame.columns:
            continue

        if target in frame.columns:
            target_series = _resolve_duplicate_column(frame, target)
            source_series = _resolve_duplicate_column(frame, source)
            frame[target] = _merge_preferred_series(target_series, source_series)
            frame = frame.drop(columns=[source])
            continue
        rename_map[source] = target
    if rename_map:
        frame = frame.rename(columns=rename_map)

    if frame.columns.duplicated().any():
        frame = frame.loc[:, ~frame.columns.duplicated()]

    for target, sources in _EXPORT_COALESCE_SOURCES.items():
        frame[target] = _coalesce_columns(frame, sources)

    for column in _EXPORT_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""

    return frame[list(_EXPORT_COLUMNS)]


def _iter_export_chunks(
    df: pd.DataFrame,
    *,
    chunk_size: int | None,
) -> Iterable[pd.DataFrame]:
    """Yield export-ready chunks with deterministic column ordering."""

    total_rows = len(df)
    if total_rows == 0:
        empty = dataframe_to_strings(df.copy(), skip=_NUMERIC_EXPORT_COLUMNS)
        yield _prepare_export_frame(empty)
        return

    effective_size = chunk_size if chunk_size and chunk_size > 0 else total_rows
    for start in range(0, total_rows, effective_size):
        stop = start + effective_size
        chunk = df.iloc[start:stop].copy()
        chunk = dataframe_to_strings(chunk, skip=_NUMERIC_EXPORT_COLUMNS)
        yield _prepare_export_frame(chunk)


def _coerce_chunk_size_value(value: object) -> int | None:
    """Return ``value`` as an ``int`` when possible, otherwise ``None``."""

    if value is None:
        return None
    if isinstance(value, bool):
        logger.warning("invalid_stream_chunk_size_bool", value=value)
        return None
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        try:
            return int(stripped)
        except ValueError:
            logger.warning("invalid_stream_chunk_size_string", value=value)
            return None
    if isinstance(value, Real):
        if pd.isna(value):
            return None
        coerced = int(value)
        if value != coerced:
            logger.warning("invalid_stream_chunk_size_float", value=float(value))
            return None
        return coerced
    if isinstance(value, pd.Series):
        if value.empty:
            return None
        if len(value) == 1:
            return _coerce_chunk_size_value(value.iloc[0])
        logger.warning(
            "invalid_stream_chunk_size_series",
            length=int(len(value)),
        )
        return None
    logger.warning(
        "invalid_stream_chunk_size_type",
        value_type=f"{type(value).__module__}.{type(value).__qualname__}",
    )
    return None


def _resolve_chunk_size(value: int | None) -> int | None:
    """Return ``value`` when positive, otherwise ``None``."""

    if value is None:
        return None
    if value <= 0:
        logger.warning("invalid_csv_chunksize", value=value)
        return None
    return value


def _resolve_stream_chunk_size(value: object) -> int:
    """Return a safe streaming chunk size derived from ``value``."""

    coerced = _coerce_chunk_size_value(value)
    resolved = _resolve_chunk_size(coerced)
    if resolved is None:
        return _EXPORT_STREAM_CHUNK_SIZE
    return resolved


def _write_export_chunks(
    chunks: Iterable[pd.DataFrame],
    path: Path,
    *,
    cfg: Config,
    key_cols: Sequence[str],
    chunk_size: int | None,
) -> Path:
    """Stream ``chunks`` to ``path`` using the deterministic CSV writer."""

    if chunk_size:
        return write_csv_chunks_deterministic(
            chunks,
            path,
            col_order=list(_EXPORT_COLUMNS),
            key_cols=list(key_cols),
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
            chunksize=chunk_size,
            merge_chunksize=chunk_size,
        )

    return write_csv_chunks_deterministic(
        chunks,
        path,
        col_order=list(_EXPORT_COLUMNS),
        key_cols=list(key_cols),
        sep=cfg.io.csv_sep,
        encoding=cfg.io.csv_encoding,
        cfg=cfg,
    )


def _maybe_run_document_postprocessing(csv_path: Path) -> None:
    if not csv_path.name.startswith("output.document_"):
        return

    data_dir: Path | None = None
    for parent in csv_path.parents:
        if parent.name.lower() == "data":
            data_dir = parent
            break

    if data_dir is None:
        return

    reference_rel = Path("input/full/document.csv")
    reference_path = data_dir / reference_rel
    if not reference_path.exists():
        return

    try:
        output_relative = csv_path.relative_to(data_dir)
    except ValueError:
        return

    ref_rel_windows = "\\".join(reference_rel.parts)
    out_rel_windows = "\\".join(output_relative.parts)
    preprocess_documents_csv(
        base_path=str(data_dir),
        ref_document_rel=ref_rel_windows,
        out_document_rel=out_rel_windows,
    )


def _finalise_export(
    df: pd.DataFrame | Iterable[pd.DataFrame],
    output: Path,
    cfg: Config,
    *,
    input_csv: Path,
    key_columns: Sequence[str] | None = None,
    chunk_size: int | None = None,
) -> int:
    """Validate input frames and write CSV/metadata artefacts."""

    if isinstance(df, pd.DataFrame):
        frames_iterable: Iterable[pd.DataFrame] = (df,)
    else:
        frames_iterable = df

    frames_iterator = iter(frames_iterable)
    required_cols = {
        name for name, col in DocumentsSchema.columns.items() if col.required
    }
    optional_cols = set(DocumentsSchema.columns) - required_cols

    present_columns: set[str] = set()
    missing_required: set[str] = set(required_cols)
    missing_optional: set[str] = set(optional_cols)

    stream_chunk = _resolve_stream_chunk_size(chunk_size)
    failure_path = output.with_name(f"{output.stem}_failure_cases.csv")
    errors = SidecarErrors()
    rows_total = 0
    rows_kept = 0
    exit_code = 0
    emitted_chunk = False
    quality_profiler = TableQualityProfiler()
    quality_summary = DocumentQualityAccumulator()

    def _validated_chunks() -> Iterator[pd.DataFrame]:
        nonlocal rows_total, rows_kept, exit_code, emitted_chunk
        nonlocal missing_required, missing_optional

        def _update_column_sets(df: pd.DataFrame) -> None:
            nonlocal missing_required, missing_optional
            present_columns.update(df.columns)
            missing_required = required_cols - present_columns
            missing_optional = optional_cols - present_columns

        for frame in frames_iterator:
            with_metadata = add_pipeline_metadata(frame)
            ordered = build_dataframe(
                with_metadata, columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False
            )
            _update_column_sets(ordered)
            rows_total += len(ordered)
            validated = ordered
            try:
                validated = DocumentsSchema.validate(
                    ordered,
                    lazy=True,
                    inplace=False,
                )
            except SchemaErrors as exc:
                for row in exc.failure_cases.to_dict("records"):
                    errors.add_error(row)
                logger.error(
                    "document_validation_failed",
                    failure_count=len(exc.failure_cases),
                    failure_path=str(failure_path),
                    error=str(exc),
                )
                validated = getattr(exc, "validated_data", ordered)
                exit_code = 1
            rows_kept += len(validated)
            cleaned = build_dataframe(
                validated, columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False
            )
            for chunk in _iter_export_chunks(cleaned, chunk_size=stream_chunk):
                emitted_chunk = True
                quality_profiler.consume(chunk)
                quality_summary.consume(chunk)
                yield chunk

        if not emitted_chunk:
            empty = build_dataframe(
                add_pipeline_metadata(pd.DataFrame()),
                columns=DOCUMENT_SCHEMA_COLUMNS,
                fill_missing=False,
            )
            _update_column_sets(empty)
            for chunk in _iter_export_chunks(empty, chunk_size=stream_chunk):
                emitted_chunk = True
                quality_profiler.consume(chunk)
                quality_summary.consume(chunk)
                yield chunk

    export_generator = _validated_chunks()

    key_cols: list[str] = []
    if key_columns:
        for column in key_columns:
            mapped = _EXPORT_COLUMN_RENAMES.get(column, column)
            if mapped in _EXPORT_COLUMNS and mapped not in key_cols:
                key_cols.append(mapped)
    if not key_cols:
        for candidate in _EXPORT_SORT_FALLBACK:
            if candidate in _EXPORT_COLUMNS:
                key_cols = [candidate]
                break
    if not key_cols:
        key_cols = [_EXPORT_COLUMNS[0]]

    col_order = list(_EXPORT_COLUMNS)
    try:
        csv_path = write_csv_chunks_deterministic(
            export_generator,
            output,
            cfg=cfg,
            key_cols=key_cols,
            col_order=col_order,
            chunksize=stream_chunk,
            merge_chunksize=stream_chunk,
            sort_chunksize=stream_chunk,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
    except OSError as exc:
        logger.error("csv_write_failed", error=str(exc), path=str(output))
        return 1

    try:
        postprocessed_path = document_export_postprocessing.postprocess_export_file(
            csv_path,
            cfg=cfg.io,
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        logger.error(
            "document_export_postprocess_failed",
            error=str(exc),
            path=str(csv_path),
        )
        exit_code = 1
    else:
        logger.info(
            "document_export_postprocess_written",
            path=str(postprocessed_path),
        )

    if missing_required:
        logger.warning(
            "validation_skipped_missing_required",
            columns=sorted(missing_required),
        )
        exit_code = 1
    elif missing_optional:
        logger.warning(
            "missing_optional_columns",
            columns=sorted(missing_optional),
        )

    errors.save(failure_path)

    rows_dropped = rows_total - rows_kept
    logger.info(
        "records_dropped",
        rows_total=int(rows_total),
        rows_kept=int(rows_kept),
        rows_dropped=int(rows_dropped),
    )
    if exit_code == 0:
        logger.info("write_done", rows=rows_kept, path=str(csv_path))
        if csv_path.name.startswith("output.document_"):
            _maybe_run_document_postprocessing(csv_path)

    doc_quality_cfg = getattr(cfg.system, "doc_quality", None)
    try:
        finalise_csv_output(
            csv_path=csv_path,
            rows_total=rows_total,
            rows_kept=rows_kept,
            command=" ".join(sys.argv),
            config_subset=_serialize_paths(cfg.to_dict()),
            inputs={"input_csv": str(input_csv)},
            schema="DocumentsSchema",
            quality_summary=quality_summary,
            quality_builder=build_quality_report,
            quality_path=csv_path.with_suffix(".quality.json"),
            quality_profiler=quality_profiler,
            quality_config=doc_quality_cfg,
            quality_table_name=csv_path.with_suffix("").name,
            quality_destination=csv_path.parent,
        )
    except QualityReportError as exc:
        destination = exc.path or csv_path.with_suffix(".quality.json")
        logger.error(
            "quality_report_write_failed",
            error=str(exc),
            path=str(destination),
        )
        return 1
    except QualityAnalysisError as exc:
        logger.exception(
            "quality_report_generation_failed",
            error=str(exc),
            exc=exc,
        )
        return 1
    return exit_code


def run_pubmed(
    cfg: Config,
    args: argparse.Namespace,
    *,
    pipeline: DocumentPipeline | None = None,
) -> int:
    """Execute the ``pubmed`` mode."""

    service = pipeline or DocumentPipeline(cfg)
    pubmed_defaults = cfg.document.pubmed
    limit = getattr(args, "limit", pubmed_defaults.limit)
    if limit is not None and limit < 0:
        logger.error(
            "invalid_limit",
            section="document.pubmed",
            limit=limit,
        )
        return 1
    offset = getattr(args, "offset", 0)
    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        legacy_output = getattr(args, "output_csv", None)
        if legacy_output not in (None, argparse.SUPPRESS):
            output_path = Path(legacy_output)
            if not isinstance(legacy_output, Path):
                args.final_out = output_path
            setattr(args, "output_csv", output_path)
        else:
            output_path = Path(io.default_output_path(args.input_csv, cfg.io))
            args.final_out = output_path
            setattr(args, "output_csv", output_path)
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        setattr(args, "output_csv", output_path)
    fallback_enabled = getattr(args, "fallback_doi_enabled", False)
    fallback_path_arg = getattr(args, "fallback_doi_path", None)
    metadata_obj = getattr(args, "_config_metadata", None)
    if not isinstance(metadata_obj, ConfigMetadata):
        metadata_obj = None
    fallback_overwrite = getattr(args, "fallback_doi_overwrite", False)
    fallback_path_text = str(fallback_path_arg) if fallback_path_arg else None
    logger.info(
        "document_pubmed_start",
        input=str(args.input_csv),
        output=str(output_path),
        limit=limit,
        offset=offset,
        workers=getattr(args, "workers", pubmed_defaults.workers),
        batch_size=getattr(args, "batch_size", pubmed_defaults.batch_size),
        sleep=getattr(args, "sleep", pubmed_defaults.sleep),
        fallback_doi_enabled=fallback_enabled,
        fallback_doi_overwrite=fallback_overwrite,
        fallback_doi_path=fallback_path_text,
        fallback_doi=prepare_option(
            metadata_obj,
            argument="fallback_doi_enabled",
            path="document.pubmed.fallback_doi_enabled",
            value=fallback_enabled,
            default_source="cli",
        ),
        fallback_doi_overwrite_meta=prepare_option(
            metadata_obj,
            argument="fallback_doi_overwrite",
            path="document.pubmed.fallback_doi_overwrite",
            value=fallback_overwrite,
            default_source="cli",
        ),
        fallback_doi_path_meta=prepare_option(
            metadata_obj,
            argument="fallback_doi_path",
            path="document.pubmed.fallback_doi_path",
            value=fallback_path_text,
            default_source="cli",
        ),
    )
    try:
        pmids_iter = io.read_ids(
            args.input_csv,
            column=getattr(args, "column", pubmed_defaults.column),
            cfg=cfg.io,
        )
    except (FileNotFoundError, ValueError) as exc:
        context = service.build_missing_input_context(Path(args.input_csv))
        logger.error(
            "input_read_failed",
            error=str(exc),
            path=str(args.input_csv),
            **context,
        )
        return 1
    if offset:
        pmids_iter = islice(pmids_iter, offset, None)
        logger.info("process_offset", offset=offset)
    pmids: Iterable[str] = pmids_iter
    limit_counter: Callable[[], int] | None = None
    if limit is not None:
        pmids_limited, get_limit_count = service.limit_iterable(pmids_iter, limit)
        pmids = pmids_limited
        limit_counter = get_limit_count

    fallback_state: FallbackDoiState | None = None
    if fallback_enabled:
        fallback_path = fallback_path_arg
        if fallback_path is None:
            logger.error(
                "fallback_doi_invalid",
                error="fallback DOI path is required when enabled",
                path="",
            )
            return 1
        delimiter = (
            getattr(args, "fallback_doi_delimiter", None) or cfg.io.csv_sep
        )
        encoding = (
            getattr(args, "fallback_doi_encoding", None) or cfg.io.csv_encoding
        )
        metrics = FallbackDoiMetrics()
        try:
            fallback_frame = pd.read_csv(
                fallback_path,
                sep=delimiter,
                encoding=encoding,
            )
        except (FileNotFoundError, pd.errors.ParserError, UnicodeError, OSError) as exc:
            logger.error(
                "fallback_doi_read_failed",
                error=str(exc),
                path=str(fallback_path),
                delimiter=delimiter,
                encoding=encoding,
            )
            return 1
        try:
            fallback_map = service.build_fallback_doi_map(
                fallback_frame,
                pmid_column=getattr(args, "fallback_doi_col_pmid", "PMID"),
                doi_column=getattr(args, "fallback_doi_col_doi", "DOI"),
                metrics=metrics,
            )
        except ValueError as exc:
            logger.error(
                "fallback_doi_invalid",
                error=str(exc),
                path=str(fallback_path),
            )
            return 1
        fallback_state = FallbackDoiState(
            path=Path(fallback_path),
            mapping=fallback_map,
            metrics=metrics,
            overwrite=getattr(args, "fallback_doi_overwrite", False),
        )

    try:
        frame_iter = service.fetch_pubmed_records(
            pmids,
            cfg=cfg,
            sleep=getattr(args, "sleep", pubmed_defaults.sleep),
            pubmed_cfg=cfg.pubmed,
            semantic_scholar_cfg=cfg.semantic_scholar,
            openalex_cfg=cfg.openalex,
            crossref_cfg=cfg.crossref,
            max_workers=getattr(args, "workers", pubmed_defaults.workers),
            batch_size=getattr(args, "batch_size", pubmed_defaults.batch_size),
            fallback_doi_map=(
                fallback_state.mapping if fallback_state else None
            ),
            return_generator=True,
        )
        output = output_path
        if fallback_state is not None:
            fallback_state.metrics.mark_total_candidates(
                len(fallback_state.mapping)
            )

            def _iter_with_fallback() -> Iterator[pd.DataFrame]:
                for frame in frame_iter:
                    yield service.apply_fallback_doi(
                        frame,
                        fallback_map=fallback_state.mapping,
                        overwrite=fallback_state.overwrite,
                        metrics=fallback_state.metrics,
                    )

            frame_iter = _iter_with_fallback()

        normalised_frames = (normalize_documents(frame) for frame in frame_iter)
        exit_code = _finalise_export(
            normalised_frames,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
            chunk_size=getattr(args, "batch_size", pubmed_defaults.batch_size),
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error(
            "pubmed_pipeline_failed",
            error=str(exc),
            output=str(output_path),
        )
        return 1
    if fallback_state is not None:
        logger.info(
            "fallback_doi_metrics",
            pipeline="pubmed",
            path=str(fallback_state.path),
            overwrite=fallback_state.overwrite,
            **fallback_state.metrics.as_log_kwargs(),
        )
    if limit_counter is not None:
        logger.info("process_limit", limit=limit_counter())
    if exit_code == 0:
        logger.info("document_pubmed_done", output=str(output_path))
    else:
        logger.error(
            "document_pubmed_failed",
            output=str(output_path),
            exit_code=exit_code,
        )
    return exit_code


def run_chembl(
    cfg: Config,
    args: argparse.Namespace,
    *,
    pipeline: DocumentPipeline | None = None,
) -> int:
    """Execute the ``chembl`` mode."""

    chembl_defaults = cfg.document.chembl
    limit = getattr(args, "limit", chembl_defaults.limit)
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="document.chembl", limit=limit)
        return 1
    offset = getattr(args, "offset", 0)
    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        legacy_output = getattr(args, "output_csv", None)
        if legacy_output not in (None, argparse.SUPPRESS):
            output_path = Path(legacy_output)
            if not isinstance(legacy_output, Path):
                args.final_out = output_path
            setattr(args, "output_csv", output_path)
        else:
            output_path = Path(io.default_output_path(args.input_csv, cfg.io))
            args.final_out = output_path
            setattr(args, "output_csv", output_path)
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        setattr(args, "output_csv", output_path)
    chunk_size = getattr(args, "chunk_size", chembl_defaults.chunk_size)
    timeout = getattr(args, "timeout", chembl_defaults.timeout)
    metadata_obj = getattr(args, "_config_metadata", None)
    if not isinstance(metadata_obj, ConfigMetadata):
        metadata_obj = None
    output_source = "cli" if getattr(args, "final_out", None) else "derived"
    service = pipeline or DocumentPipeline(cfg)
    logger.info(
        "document_chembl_start",
        input=prepare_option(metadata_obj, value=str(args.input_csv), default_source="cli"),
        output=prepare_option(
            metadata_obj,
            value=str(output_path),
            default_source=output_source,
        ),
        limit=prepare_option(
            metadata_obj,
            argument="limit",
            path="sources.chembl.pipelines.document.chembl.limit",
            value=limit,
        ),
        offset=prepare_option(
            metadata_obj,
            argument="offset",
            path="sources.chembl.pipelines.document.chembl.offset",
            value=offset,
            default_source="cli",
        ),
        chunk_size=prepare_option(
            metadata_obj,
            argument="chunk_size",
            path="sources.chembl.pipelines.document.chembl.chunk_size",
            value=chunk_size,
        ),
        timeout=prepare_option(
            metadata_obj,
            argument="timeout",
            path="sources.chembl.pipelines.document.chembl.timeout",
            value=timeout,
        ),
    )

    rate_cfg = cfg.rate
    global_limiter = None
    if (rate_cfg.global_rps or 0) > 0:
        global_limiter = get_global_limiter(rate_cfg.global_rps, rate_cfg.global_burst)

    with ChemblClient(
        cfg.api, cfg.retry, cfg.chembl, global_limiter=global_limiter
    ) as client:
        try:
            ids_iter = io.read_ids(
                args.input_csv,
                column=getattr(args, "column", chembl_defaults.column),
                cfg=cfg.io,
            )
        except (FileNotFoundError, ValueError) as exc:
            context = service.build_missing_input_context(Path(args.input_csv))
            logger.error(
                "input_read_failed",
                error=str(exc),
                path=str(args.input_csv),
                **context,
            )
            return 1

        if offset:
            ids_iter = islice(ids_iter, offset, None)
            logger.info("process_offset", offset=offset)

        ids: Iterable[str] = ids_iter
        limit_counter: Callable[[], int] | None = None
        if limit is not None:
            limited_ids, get_limit_count = service.limit_iterable(ids_iter, limit)
            ids = limited_ids
            limit_counter = get_limit_count

        try:
            df = cl.get_documents(
                ids,
                cfg=cfg.api,
                client=client,
                chunk_size=getattr(args, "chunk_size", chembl_defaults.chunk_size),
                timeout=getattr(args, "timeout", chembl_defaults.timeout),
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error(
                "chembl_documents_fetch_failed",
                error=str(exc),
                chunk_size=getattr(args, "chunk_size", chembl_defaults.chunk_size),
            )
            return 1
        if "doi" in df.columns:
            df["doi"] = df["doi"].map(normalise_doi)
        output = output_path
        df = normalize_documents(df)
        exit_code = _finalise_export(
            df,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
            chunk_size=getattr(args, "chunk_size", chembl_defaults.chunk_size),
        )
        if exit_code == 0:
            logger.info("document_chembl_done", output=str(output_path))
        else:
            logger.error(
                "document_chembl_failed",
                output=str(output_path),
                exit_code=exit_code,
            )
        if limit_counter is not None:
            logger.info("process_limit", limit=limit_counter())
        return exit_code


def run_all(
    cfg: Config,
    args: argparse.Namespace,
    *,
    pipeline: DocumentPipeline | None = None,
) -> int:
    """Execute the ``all`` mode by merging ChEMBL and PubMed outputs."""

    service = pipeline or DocumentPipeline(cfg)
    all_defaults = cfg.document.all
    limit = getattr(args, "limit", all_defaults.limit)
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="document.all", limit=limit)
        return 1

    rate_cfg = cfg.rate
    global_limiter = None
    if (rate_cfg.global_rps or 0) > 0:
        global_limiter = get_global_limiter(rate_cfg.global_rps, rate_cfg.global_burst)

    try:
        ids_iter = io.read_ids(
            args.input_csv,
            column=getattr(args, "column", all_defaults.column),
            cfg=cfg.io,
        )
    except (FileNotFoundError, ValueError) as exc:
        context = service.build_missing_input_context(Path(args.input_csv))
        logger.error(
            "input_read_failed",
            error=str(exc),
            path=str(args.input_csv),
            **context,
        )
        return 1

    offset = getattr(args, "offset", 0)
    if offset:
        ids_iter = islice(ids_iter, offset, None)
        logger.info("process_offset", offset=offset)

    ids_source: Iterable[str] = ids_iter
    limit_counter: Callable[[], int] | None = None
    if limit is not None:
        ids_limited, get_limit_count = service.limit_iterable(ids_iter, limit)
        ids_source = ids_limited
        limit_counter = get_limit_count

    iterator = iter(ids_source)
    sample_size = getattr(args, "chembl_chunk_size", all_defaults.chunk_size)
    sample_ids = list(islice(iterator, sample_size))
    ids_for_fetch = chain(sample_ids, iterator)
    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        legacy_output = getattr(args, "output_csv", None)
        if legacy_output not in (None, argparse.SUPPRESS):
            output_path = Path(legacy_output)
            if not isinstance(legacy_output, Path):
                args.final_out = output_path
            setattr(args, "output_csv", output_path)
        else:
            output_path = Path(io.default_output_path(args.input_csv, cfg.io))
            args.final_out = output_path
            setattr(args, "output_csv", output_path)
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        setattr(args, "output_csv", output_path)
    fallback_enabled = getattr(args, "fallback_doi_enabled", False)
    fallback_path_arg = getattr(args, "fallback_doi_path", None)
    logger.info(
        "document_all_start",
        input=str(args.input_csv),
        output=str(output_path),
        limit=limit,
        offset=offset,
        chembl_chunk_size=getattr(args, "chembl_chunk_size", all_defaults.chunk_size),
        pubmed_workers=getattr(args, "pubmed_workers", all_defaults.workers),
        pubmed_batch_size=getattr(args, "pubmed_batch_size", all_defaults.batch_size),
        pubmed_sleep=getattr(args, "pubmed_sleep", all_defaults.sleep),
        chembl_timeout=getattr(args, "chembl_timeout", all_defaults.timeout),
        fallback_doi_enabled=fallback_enabled,
        fallback_doi_overwrite=getattr(args, "fallback_doi_overwrite", False),
        fallback_doi_path=str(fallback_path_arg) if fallback_path_arg else None,
    )

    fallback_state: FallbackDoiState | None = None
    if fallback_enabled:
        fallback_path = fallback_path_arg
        if fallback_path is None:
            logger.error(
                "fallback_doi_invalid",
                error="fallback DOI path is required when enabled",
                path="",
            )
            return 1
        delimiter = (
            getattr(args, "fallback_doi_delimiter", None) or cfg.io.csv_sep
        )
        encoding = (
            getattr(args, "fallback_doi_encoding", None) or cfg.io.csv_encoding
        )
        metrics = FallbackDoiMetrics()
        try:
            fallback_frame = pd.read_csv(
                fallback_path,
                sep=delimiter,
                encoding=encoding,
            )
        except (FileNotFoundError, pd.errors.ParserError, UnicodeError, OSError) as exc:
            logger.error(
                "fallback_doi_read_failed",
                error=str(exc),
                path=str(fallback_path),
                delimiter=delimiter,
                encoding=encoding,
            )
            return 1
        try:
            fallback_map = service.build_fallback_doi_map(
                fallback_frame,
                pmid_column=getattr(args, "fallback_doi_col_pmid", "PMID"),
                doi_column=getattr(args, "fallback_doi_col_doi", "DOI"),
                metrics=metrics,
            )
        except ValueError as exc:
            logger.error(
                "fallback_doi_invalid",
                error=str(exc),
                path=str(fallback_path),
            )
            return 1
        fallback_state = FallbackDoiState(
            path=Path(fallback_path),
            mapping=fallback_map,
            metrics=metrics,
            overwrite=getattr(args, "fallback_doi_overwrite", False),
        )
        fallback_state.metrics.mark_total_candidates(len(fallback_state.mapping))

    try:
        with ChemblClient(
            cfg.api, cfg.retry, cfg.chembl, global_limiter=global_limiter
        ) as client:
            doc_df = cl.get_documents(
                ids_for_fetch,
                cfg=cfg.api,
                client=client,
                chunk_size=getattr(
                    args, "chembl_chunk_size", all_defaults.chunk_size
                ),
                timeout=getattr(
                    args, "chembl_timeout", all_defaults.timeout
                ),
            )
    except (requests.RequestException, ValueError) as exc:
        logger.error(
            "chembl_documents_fetch_failed",
            ids=sample_ids,
            error=str(exc),
            output=str(output_path),
            chunk_size=getattr(args, "chembl_chunk_size", all_defaults.chunk_size),
        )
        return 1
    if limit_counter is not None:
        logger.info("process_limit", limit=limit_counter())
    output = output_path
    if "doi" in doc_df.columns:
        doc_df["doi"] = doc_df["doi"].map(normalise_doi)
    if fallback_state is not None:
        doc_df = service.apply_fallback_doi(
            doc_df,
            fallback_map=fallback_state.mapping,
            overwrite=fallback_state.overwrite,
            metrics=fallback_state.metrics,
            pmid_column="pubmed_id",
        )
    if doc_df.empty or "pubmed_id" not in doc_df:
        processed = dp.postprocess_documents(doc_df)
        extra_cols = [c for c in doc_df.columns if c not in processed.columns]
        if extra_cols:
            processed = processed.merge(
                doc_df[["document_chembl_id"] + extra_cols],
                on="document_chembl_id",
                how="left",
            )
        processed = normalize_documents(processed)
        exit_code = _finalise_export(
            processed,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
            chunk_size=getattr(
                args, "chembl_chunk_size", all_defaults.chunk_size
            ),
        )
        if fallback_state is not None:
            logger.info(
                "fallback_doi_metrics",
                pipeline="all",
                path=str(fallback_state.path),
                overwrite=fallback_state.overwrite,
                **fallback_state.metrics.as_log_kwargs(),
            )
        if exit_code == 0:
            logger.info("document_all_done", output=str(output_path))
        else:
            logger.error(
                "document_all_failed",
                output=str(output_path),
                exit_code=exit_code,
            )
        return exit_code

    pubmed_ids = pd.to_numeric(doc_df["pubmed_id"], errors="coerce").astype("Int64")
    doi_fallback_map: dict[str, str] = {}
    if "doi" in doc_df.columns:
        doi_series = doc_df["doi"].astype("string")
        mask = pubmed_ids.notna() & doi_series.notna() & (doi_series != "")
        if mask.any():
            masked_pmids = pubmed_ids[mask].tolist()
            masked_dois = doi_series[mask].tolist()
            if len(masked_pmids) != len(masked_dois):
                raise ValueError("mismatched DOI fallback map lengths")
            doi_fallback_map = {
                str(pmid): normalise_doi(doi)
                for pmid, doi in zip(masked_pmids, masked_dois, strict=True)
            }
    if fallback_state is not None and fallback_state.mapping:
        combined_map = dict(doi_fallback_map)
        combined_map.update(fallback_state.mapping)
        doi_fallback_map = combined_map
    pmids = pubmed_ids.dropna().astype(str).tolist()
    pubmed_batch_size = getattr(args, "pubmed_batch_size", all_defaults.batch_size)
    if pubmed_batch_size is None or pubmed_batch_size <= 0:
        pubmed_batch_size = all_defaults.batch_size
    merge_chunk_size = getattr(args, "chembl_chunk_size", all_defaults.chunk_size)
    if merge_chunk_size is None or merge_chunk_size <= 0:
        merge_chunk_size = all_defaults.chunk_size

    pubmed_frames = service.fetch_pubmed_records(
        pmids,
        cfg=cfg,
        sleep=getattr(args, "pubmed_sleep", all_defaults.sleep),
        semantic_scholar_cfg=cfg.semantic_scholar,
        openalex_cfg=cfg.openalex,
        crossref_cfg=cfg.crossref,
        max_workers=getattr(args, "pubmed_workers", all_defaults.workers),
        batch_size=pubmed_batch_size,
        pubmed_cfg=cfg.pubmed,
        fallback_doi_map=doi_fallback_map or None,
        return_generator=True,
    )
    doc_df["pubmed_id"] = pubmed_ids.astype("Int64").astype("string").fillna("")

    try:
        with tempfile.TemporaryDirectory(prefix="chembl_pubmed_") as tmp_dir:
            tmp_path = Path(tmp_dir) / "pubmed_metadata.csv"
            metadata_path = write_csv_chunks_deterministic(
                pubmed_frames,
                tmp_path,
                key_cols=["PubMed.PMID"],
                chunksize=pubmed_batch_size,
                merge_chunksize=pubmed_batch_size,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                cfg=cfg,
            )
            if metadata_path.exists() and metadata_path.stat().st_size > 0:
                metadata_iter = service.read_csv_chunks(
                    metadata_path,
                    cfg=cfg,
                    chunk_size=merge_chunk_size,
                )
            else:
                metadata_iter = iter(())
            merged = merge_with_chembl(doc_df, metadata_iter)
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error(
            "pubmed_pipeline_failed",
            error=str(exc),
            output=str(output_path),
        )
        return 1
    processed = dp.postprocess_documents(merged)
    extra_cols = [c for c in merged.columns if c not in processed.columns]
    if extra_cols:
        processed = processed.merge(
            merged[["document_chembl_id"] + extra_cols],
            on="document_chembl_id",
            how="left",
        )
    processed = normalize_documents(processed)
    exit_code = _finalise_export(
        processed,
        output,
        cfg,
        input_csv=Path(args.input_csv),
        key_columns=["document_chembl_id"],
        chunk_size=getattr(args, "chembl_chunk_size", all_defaults.chunk_size),
    )
    if fallback_state is not None:
        logger.info(
            "fallback_doi_metrics",
            pipeline="all",
            path=str(fallback_state.path),
            overwrite=fallback_state.overwrite,
            **fallback_state.metrics.as_log_kwargs(),
        )
    if exit_code == 0:
        logger.info("document_all_done", output=str(output_path))
    else:
        logger.error(
            "document_all_failed",
            output=str(output_path),
            exit_code=exit_code,
        )
    return exit_code


MODE_HANDLERS: Mapping[
    str, Callable[[Config, argparse.Namespace, DocumentPipeline], int]
] = {
    "chembl": run_chembl,
    "pubmed": run_pubmed,
    "all": run_all,
}


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the selected document pipeline with CLI-specific hooks."""

    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        legacy_output = getattr(args, "output_csv", None)
        if legacy_output not in (None, argparse.SUPPRESS):
            output_path = Path(legacy_output)
            if not isinstance(legacy_output, Path):
                args.final_out = output_path
        else:
            output_path = Path(io.default_output_path(args.input_csv, cfg.io))
            args.final_out = output_path
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
    setattr(args, "output_csv", output_path)
    mode = getattr(args, "mode", None)
    if mode in (None, ""):
        mode = getattr(args, "command", None)
    handler = MODE_HANDLERS.get(str(mode))
    if handler is None:
        logger.error(
            "missing_subcommand_handler",
            command=str(mode) if mode is not None else "",
        )
        return 1
    if getattr(args, "skip_existing", False) and not getattr(args, "force", False):
        output_path = Path(getattr(args, "final_out", output_path))
        if output_path.exists():
            logger.info("pipeline_skip_existing", output=str(output_path))
            return 0
    timeout = getattr(args, "timeout", None)
    if timeout is not None:
        cfg.api.timeout_read = timeout
        cfg.api.timeout_connect = timeout
    service = DocumentPipeline(cfg)
    current_handler = handler

    supports_pipeline = True
    try:
        signature = inspect.signature(current_handler)
    except (TypeError, ValueError):
        signature = None

    if signature is not None:
        supports_pipeline = any(
            parameter.kind is inspect.Parameter.VAR_KEYWORD
            or parameter.name == "pipeline"
            for parameter in signature.parameters.values()
            if parameter.kind
            in {
                inspect.Parameter.POSITIONAL_OR_KEYWORD,
                inspect.Parameter.KEYWORD_ONLY,
                inspect.Parameter.VAR_KEYWORD,
            }
        )

    if supports_pipeline:
        result = current_handler(cfg, args, pipeline=service)
    else:
        result = current_handler(cfg, args)
    return int(result)


__all__ = [
    "MODE_HANDLERS",
    "run",
    "run_all",
    "run_chembl",
    "run_pubmed",
]
