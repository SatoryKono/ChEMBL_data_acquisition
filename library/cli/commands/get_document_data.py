"""Command line interface for retrieving document metadata from external sources.

The tool integrates :mod:`library.integration.pubmed_library` and
:mod:`library.integration.chembl_library` to collect information about publications from
several public APIs.  The interface mirrors :mod:`scripts.get_target_data` and exposes a
single entry point configured via ``--mode``:

``--mode pubmed``
    Query PubMed, Semantic Scholar, OpenAlex and CrossRef for a list of PMIDs.
``--mode chembl``
    Retrieve document information from the ChEMBL API.
``--mode all``
    Run the ChEMBL and PubMed pipelines and merge the results.

Example
-------
Fetch PubMed metadata for identifiers listed in ``pmids.csv``::

    python scripts/get_document_data.py --mode pubmed --config config/config.yaml --input pmids.csv --final-out output.csv

The input file must contain a ``PMID`` column.

By default the command writes the canonical dataset together with quality-control
and correlation summaries using :func:`library.io.save_standard_outputs`.  The
artefacts are timestamped with the UTC date in ``YYYYMMDD`` format and stored in
``cfg.io.output_dir``.  Legacy sidecar files (failure CSVs, metadata YAML and
``.quality.json``) are only emitted when ``--emit-legacy-artifacts`` is
enabled.

"""

from __future__ import annotations

import argparse
import inspect
import logging
import os
import re
import sys
import tempfile
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from itertools import chain, islice
from numbers import Integral, Real
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from pandera.errors import SchemaErrors

from library import io
from library.cli import (
    ConfigMetadata,
    Logger,
    LoggerConfig,
    build_root_parser,
    configure_logger,
    path_argument,
    positive_int,
    prepare_io_paths,
    set_emit_legacy_help,
)
from library.cli.logging import setup_cli_logging
from library.cli.metadata import prepare_option
from library.cli.utils import run_cli_command
from library.common.csv_utils import write_csv_chunks_deterministic
from library.common.log import logger
from library.common.run_context import get_current as get_run_context
from library.common.sidecar import SidecarErrors, resolve_failure_chunk_size
from library.config import (
    Config,
    _serialize_paths,
)
from library.document_defaults import ALL_DEFAULTS, CHEMBL_DEFAULTS, PUBMED_DEFAULTS
from library.integration import chembl_library as cl
from library.orchestration import ETLContext
from library.pipelines.common import add_pipeline_metadata
from library.pipelines.common.metadata import get_pipeline_version
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
from library.qa.reporting import build_table_quality_hook
from library.qa.table_quality import TableQualityProfiler
from library.reporting.run_manifest import (
    QualityAnalysisError,
    QualityReportError,
    finalise_csv_output,
)
from library.schemas import DocumentsSchema, normalize_documents
from library.schemas.document_spec import DOCUMENT_EXPORT_COLUMNS
from library.utils.data_correlation import generate_correlation_report
from library.utils.qc_report import generate_qc_report
from library.validation import validate_documents

DEFAULT_INPUT_NAME = "document.csv"
DEFAULT_OUTPUT_STEM = "documents"

_CHEMBL_METADATA_SOURCES: tuple[str, ...] = ("ChEMBL Document API",)
_PUBMED_METADATA_SOURCES: tuple[str, ...] = (
    "PubMed",
    "Semantic Scholar",
    "OpenAlex",
    "CrossRef",
)
_ALL_METADATA_SOURCES: tuple[str, ...] = (
    *_CHEMBL_METADATA_SOURCES,
    *_PUBMED_METADATA_SOURCES,
)


class _FallbackPathAction(argparse.Action):
    """Store fallback DOI CSV path under both legacy and new attribute names."""

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: object,
        option_string: str | None = None,
    ) -> None:
        setattr(namespace, self.dest, values)
        namespace.fallback_doi_csv = values


# ``DOCUMENT_EXPORT_COLUMNS`` is provided as an immutable tuple in the schema
# declaration to make accidental mutations unlikely.  Pandas, however, expects a
# list-like object for column projections and interprets tuples as single column
# keys (raising ``KeyError`` on recent releases).  Store the export projection as
# a list locally to keep the deterministic ordering while remaining compatible
# with pandas' expectations.
_EXPORT_COLUMNS: list[str] = list(DOCUMENT_EXPORT_COLUMNS)


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

_TABLE_NAME_PREFIX = "output."
_DATE_SUFFIX_RE = re.compile(r"_(\d{8})$")


def _resolve_table_name_and_date(output: Path) -> tuple[str, str | None]:
    """Return a normalised table name and optional date inferred from ``output``."""

    stem = output.stem or DEFAULT_OUTPUT_STEM
    inferred_date: str | None = None
    match = _DATE_SUFFIX_RE.search(stem)
    if match is not None:
        candidate = match.group(1)
        if candidate.isdigit():
            inferred_date = candidate
            stem = stem[: match.start()]

    if stem.startswith(_TABLE_NAME_PREFIX) and len(stem) > len(_TABLE_NAME_PREFIX):
        stem = stem[len(_TABLE_NAME_PREFIX) :]

    normalised = stem.strip("._") or DEFAULT_OUTPUT_STEM
    return normalised, inferred_date


def _resolve_timeout(value: float | None, default: float) -> float:
    """Return ``default`` when ``value`` is ``None`` otherwise the float value."""

    if value is None:
        return float(default)
    return float(value)


def _iter_export_chunks(df: pd.DataFrame, *, chunk_size: int) -> Iterable[pd.DataFrame]:
    """Yield export-ready DataFrame chunks from ``df``."""

    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    if df.empty:
        yield build_dataframe([], columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False)
        return

    total = len(df)
    for start in range(0, total, chunk_size):
        stop = start + chunk_size
        chunk = df.iloc[start:stop]
        export_chunk = build_dataframe(
            chunk, columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False
        )
        export_chunk = dataframe_to_strings(export_chunk, skip=_NUMERIC_EXPORT_COLUMNS)
        yield _prepare_export_frame(export_chunk)


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
        selected.replace("", pd.NA).bfill(axis=1).ffill(axis=1).iloc[:, 0].fillna("")
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


def _merge_preferred_series(
    target_series: pd.Series, source_series: pd.Series
) -> pd.Series:
    """Return ``target_series`` with missing entries populated from ``source_series``."""

    if target_series.empty and source_series.empty:
        return target_series.copy()

    if not target_series.index.equals(source_series.index):
        source_series = source_series.reindex(target_series.index)

    combined = target_series.copy()

    missing_mask = target_series.isna()
    if missing_mask.any():
        combined.loc[missing_mask] = source_series.loc[missing_mask]

    if pd.api.types.is_object_dtype(
        target_series.dtype
    ) or pd.api.types.is_string_dtype(target_series.dtype):
        empty_mask = target_series.fillna("").eq("")
        if empty_mask.any():
            combined.loc[empty_mask] = source_series.loc[empty_mask]

    return combined


def _prepare_export_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Rename and project columns to match the export schema."""

    # Coalesce legacy column names into the canonical ``ChEMBL.*`` aliases while
    # keeping existing data intact.

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

    # ``DataFrame.reindex`` guarantees that all expected columns exist while keeping the
    # deterministic ordering defined by ``_EXPORT_COLUMNS``.  Missing columns are filled
    # with empty strings which mirrors the legacy behaviour of manually initialising
    # them in a loop.  Returning the reindexed frame avoids the ``KeyError`` raised by
    # ``frame[_EXPORT_COLUMNS]`` on newer pandas releases where a tuple of column names
    # is interpreted as a single label rather than a sequence.
    return frame.reindex(columns=list(_EXPORT_COLUMNS), fill_value="")


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


@dataclass(slots=True)
class FinaliseExportResult:
    """Outcome of ``_finalise_export`` including postprocess artefact details."""

    exit_code: int
    artifacts: io.StandardOutputArtifacts | None = None
    legacy_csv_path: Path | None = None
    failure_path: Path | None = None
    postprocess_path: Path | None = None
    postprocess_metrics: object | None = None
    postprocess_report: Path | None = None
    date_tag: str | None = None


def _run_documents_postprocess(
    csv_path: Path,
    destination: Path,
    *,
    config_override: Path | str | None = None,
    extras: Mapping[str, object] | None = None,
) -> tuple[Path, object | None, Path | None]:
    """Execute the documents postprocess pipeline and return artefact details."""

    from library.postprocess.common import (
        PostprocessingPipelineConfig,
        generate_metrics_report,
        get_csv_runtime_config,
        get_pipeline_config,
        run_postprocessing_pipeline,
    )
    from library.postprocessing.documents import (
        DOCUMENT_SCHEMA as POSTPROCESS_DOCUMENT_SCHEMA,
    )
    from library.postprocessing.documents import (
        run_document_pipeline as run_document_postprocess,
    )
    from library.postprocessing.documents import steps as document_steps
    from library.postprocessing.documents import (
        validate_documents as validate_document_postprocess,
    )

    pipeline_config = get_pipeline_config("documents", config_override)
    csv_cfg = get_csv_runtime_config(pipeline_config)

    document_steps.PIPELINE_CONFIG = pipeline_config
    document_steps.PIPELINE_STEPS = pipeline_config.step_definitions()

    runtime = PostprocessingPipelineConfig(
        pipeline_config=pipeline_config,
        csv_runtime_config=csv_cfg,
        runner=run_document_postprocess,
        validator=validate_document_postprocess,
        schema=POSTPROCESS_DOCUMENT_SCHEMA,
        logger=logger,
    )

    result = run_postprocessing_pipeline(
        "documents",
        csv_path,
        destination,
        runtime,
    )

    metrics = result.metrics
    report_path = result.report_path
    pipeline_version = (
        metrics.pipeline_version if metrics and metrics.pipeline_version else None
    ) or pipeline_config.pipeline_version

    if extras or metrics is None or report_path is None:
        metrics, report_path = generate_metrics_report(
            "documents",
            destination,
            csv_cfg,
            run_document_postprocess,
            pipeline_version=pipeline_version,
            extras=extras,
            logger=logger,
            pipeline_metrics=metrics,
        )

    return result.output_path, metrics, report_path


def _finalise_export(
    df: pd.DataFrame | Iterable[pd.DataFrame],
    output: Path,
    cfg: Config,
    *,
    input_csv: Path,
    key_columns: Sequence[str] | None = None,
    chunk_size: int | None = None,
    partial_run: bool = False,
    rerun_postprocess: bool = False,
    postprocess_enabled: bool = False,
    postprocess_config_path: Path | str | None = None,
    emit_legacy_artifacts: bool = False,
    date_tag: str | None = None,
    cli_args: argparse.Namespace | None = None,
    metadata_sources: Sequence[str] | None = None,
    stats_extra: Mapping[str, Any] | None = None,
) -> FinaliseExportResult:
    """Validate input frames and persist standard pipeline artefacts."""

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
    errors = SidecarErrors(chunk_size=resolve_failure_chunk_size(cfg))
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
                validation = validate_documents(ordered, return_result=True)
            except SchemaErrors as exc:
                for row in exc.failure_cases.to_dict("records"):
                    errors.add_error(row)
                logger.error(
                    "document_validation_failed",
                    failure_count=len(exc.failure_cases),
                    failure_path=str(failure_path),
                    error=str(exc),
                    exc_info=exc,
                )
                validated = getattr(exc, "validated_data", ordered)
                exit_code = 1
            else:
                validated = validation.data
                if not validation.failure_cases.empty:
                    failure_records = validation.failure_cases.to_dict("records")
                    for row in failure_records:
                        errors.add_error(row)
                    logger.error(
                        "document_validation_failed",
                        failure_count=len(validation.failure_cases),
                        failure_path=str(failure_path),
                    )
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

    validated_chunks = list(_validated_chunks())

    if validated_chunks:
        if len(validated_chunks) == 1:
            export_frame = validated_chunks[0].copy()
        else:
            export_frame = pd.concat(validated_chunks, ignore_index=True)
            export_frame = export_frame.reindex(
                columns=list(_EXPORT_COLUMNS), fill_value=""
            )
    else:  # pragma: no cover - defensive guard
        export_frame = build_dataframe(
            add_pipeline_metadata(pd.DataFrame()),
            columns=DOCUMENT_SCHEMA_COLUMNS,
            fill_missing=False,
        )
        export_frame = dataframe_to_strings(export_frame, skip=_NUMERIC_EXPORT_COLUMNS)
        export_frame = _prepare_export_frame(export_frame)

    table_name, inferred_date = _resolve_table_name_and_date(output)
    resolved_date_tag = (
        date_tag or inferred_date or datetime.now(UTC).strftime("%Y%m%d")
    )

    try:
        # Добавляем индексную колонку raw.index в основную таблицу
        if "raw.index" not in export_frame.columns:
            export_frame.insert(0, "raw.index", export_frame.index)
            logging.info(
                "Добавлена индексная колонка 'raw.index' (%s строк).",
                len(export_frame),
            )

        quality_report = generate_qc_report(
            export_frame,
            table_name=table_name,
            profiler=quality_profiler,
        )
        correlation_report = generate_correlation_report(
            export_frame,
            table_name=table_name,
            profiler=quality_profiler,
        )
    except ValueError as exc:  # pragma: no cover - validation guard
        logger.error(
            "document_quality_build_failed",
            error=str(exc),
            exc_info=exc,
            table=table_name,
        )
        return FinaliseExportResult(exit_code=1, date_tag=resolved_date_tag)

    try:
        artifacts = io.save_standard_outputs(
            export_frame,
            correlation_report,
            quality_report,
            table_name=table_name,
            date_tag=resolved_date_tag,
            output_path=output,
        )
    except (OSError, ValueError) as exc:
        logger.error(
            "document_standard_outputs_failed",
            error=str(exc),
            exc_info=exc,
            table=table_name,
        )
        return FinaliseExportResult(exit_code=1, date_tag=resolved_date_tag)

    csv_path = artifacts.dataset

    qc_summary_payload = quality_summary.build()
    stats_payload = dict(stats_extra) if stats_extra else None

    io.save_metadata(
        table_name=table_name,
        date_tag=resolved_date_tag,
        args=cli_args,
        qc_summary=qc_summary_payload,
        output_dir=csv_path.parent,
        artifacts=[
            artifacts.dataset,
            artifacts.quality_report,
            artifacts.correlation_report,
        ],
        sources=metadata_sources,
        run_context=get_run_context(),
        stats_extra=stats_payload,
    )

    logger.info(
        "document_standard_outputs_written",
        dataset=str(artifacts.dataset),
        quality_report=str(artifacts.quality_report),
        correlation_report=str(artifacts.correlation_report),
        date_tag=resolved_date_tag,
    )

    postprocessed_path: Path | None = None
    postprocess_metrics: object | None = None
    postprocess_report: Path | None = None

    if not postprocess_enabled:
        logger.info("[INFO] Postprocessing skipped (flag --postprocess not set)")
    else:
        postprocess_extras: Mapping[str, object] | None = None
        if emit_legacy_artifacts and partial_run:
            postprocess_extras = {"partial_run": True}

        destination = csv_path.with_name("output_postprocessed.documents.csv")

        try:
            from library.postprocessing.common.types import (
                SchemaValidationError,
                StepError,
            )
        except ImportError:  # pragma: no cover - defensive guard
            SchemaValidationError = StepError = RuntimeError  # type: ignore[assignment]

        try:
            postprocessed_path, postprocess_metrics, postprocess_report = (
                _run_documents_postprocess(
                    csv_path,
                    destination,
                    config_override=postprocess_config_path,
                    extras=postprocess_extras,
                )
            )
        except (  # type: ignore[misc]
            SchemaValidationError,
            StepError,
            OSError,
            ValueError,
            pd.errors.ParserError,
        ) as exc:
            logger.error(
                "document_export_postprocess_failed",
                error=str(exc),
                exc_info=exc,
                path=str(csv_path),
            )
            postprocessed_path = None
            postprocess_metrics = None
            postprocess_report = None
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

    legacy_failure_path: Path | None = None
    legacy_csv_path: Path | None = None

    if emit_legacy_artifacts:
        errors.save(failure_path, cfg=cfg)
        legacy_failure_path = failure_path if failure_path.exists() else None

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

        try:
            legacy_csv_path = _write_export_chunks(
                validated_chunks,
                output,
                cfg=cfg,
                key_cols=key_cols,
                chunk_size=stream_chunk,
            )
        except OSError as exc:
            logger.error(
                "csv_write_failed",
                error=str(exc),
                exc_info=exc,
                path=str(output),
            )
            exit_code = 1
            legacy_csv_path = None
        else:
            logger.info("legacy_write_done", rows=rows_kept, path=str(legacy_csv_path))

        if legacy_csv_path is not None:
            doc_quality_cfg = getattr(cfg.system, "doc_quality", None)
            quality_hook = build_table_quality_hook(
                doc_quality_cfg,
                table_name=legacy_csv_path.with_suffix(""),
                destination=legacy_csv_path.parent,
            )
            try:
                finalise_csv_output(
                    csv_path=legacy_csv_path,
                    rows_total=rows_total,
                    rows_kept=rows_kept,
                    command=" ".join(sys.argv),
                    config=_serialize_paths(cfg.to_dict()),
                    inputs={"input_csv": str(input_csv)},
                    schema="DocumentsSchema",
                    stats_extra=stats_payload,
                    quality_summary=quality_summary,
                    quality_builder=build_quality_report,
                    quality_path=legacy_csv_path.with_suffix(".quality.json"),
                    quality_profiler=quality_profiler,
                    quality_hook=quality_hook,
                )
            except QualityReportError as exc:
                destination = exc.path or legacy_csv_path.with_suffix(".quality.json")
                logger.error(
                    "quality_report_write_failed",
                    error=str(exc),
                    exc_info=exc,
                    path=str(destination),
                )
                return FinaliseExportResult(
                    exit_code=1,
                    artifacts=artifacts,
                    legacy_csv_path=legacy_csv_path,
                    failure_path=legacy_failure_path,
                    postprocess_path=postprocessed_path,
                    postprocess_metrics=postprocess_metrics,
                    postprocess_report=postprocess_report,
                    date_tag=resolved_date_tag,
                )
            except QualityAnalysisError as exc:
                logger.exception(
                    "quality_report_generation_failed",
                    error=str(exc),
                    exc=exc,
                )
                return FinaliseExportResult(
                    exit_code=1,
                    artifacts=artifacts,
                    legacy_csv_path=legacy_csv_path,
                    failure_path=legacy_failure_path,
                    postprocess_path=postprocessed_path,
                    postprocess_metrics=postprocess_metrics,
                    postprocess_report=postprocess_report,
                    date_tag=resolved_date_tag,
                )
    else:
        logger.debug(
            "legacy_artifacts_skipped",
            flag="--emit-legacy-artifacts",
            output=str(output),
        )

    rows_dropped = rows_total - rows_kept
    logger.info(
        "records_dropped",
        rows_total=int(rows_total),
        rows_kept=int(rows_kept),
        rows_dropped=int(rows_dropped),
    )
    if exit_code == 0:
        logger.info("write_done", rows=rows_kept, path=str(csv_path))

    return FinaliseExportResult(
        exit_code=exit_code,
        artifacts=artifacts,
        legacy_csv_path=legacy_csv_path,
        failure_path=legacy_failure_path,
        postprocess_path=postprocessed_path,
        postprocess_metrics=postprocess_metrics,
        postprocess_report=postprocess_report,
        date_tag=resolved_date_tag,
    )


def _log_document_completion(
    event: str,
    *,
    cfg: Config,
    output_path: Path,
    logger: Logger,
    extras: Mapping[str, object] | None = None,
    finalise_result: FinaliseExportResult | None = None,
    postprocess_enabled: bool = False,
) -> None:
    """Log pipeline completion details and write the postprocess report."""

    extras_payload: dict[str, object] = dict(extras) if extras else {}
    metrics = None
    postprocess_path: Path | None = None
    report_path: Path | None = None
    dataset_path: Path | None = None
    quality_report_path: Path | None = None
    correlation_path: Path | None = None
    legacy_csv_path: Path | None = None
    date_tag = None

    if finalise_result is not None:
        postprocess_path = finalise_result.postprocess_path
        metrics = finalise_result.postprocess_metrics
        report_path = finalise_result.postprocess_report
        legacy_csv_path = finalise_result.legacy_csv_path
        date_tag = finalise_result.date_tag
        if finalise_result.artifacts is not None:
            dataset_path = finalise_result.artifacts.dataset
            quality_report_path = finalise_result.artifacts.quality_report
            correlation_path = finalise_result.artifacts.correlation_report

    if postprocess_enabled and postprocess_path is not None:
        extras_payload["postprocess_output"] = str(postprocess_path)
    if postprocess_enabled and report_path is not None:
        extras_payload["postprocess_report"] = str(report_path)

    pipeline_version_value = get_pipeline_version()
    if metrics is not None and getattr(metrics, "pipeline_version", None) is not None:
        pipeline_version_value = metrics.pipeline_version

    primary_output = postprocess_path or dataset_path or output_path

    payload: dict[str, object] = {
        "output_postprocessed": str(primary_output),
        "pipeline_version": pipeline_version_value,
    }
    if dataset_path is not None:
        payload["standard_dataset"] = str(dataset_path)
    if quality_report_path is not None:
        payload["standard_quality_report"] = str(quality_report_path)
    if correlation_path is not None:
        payload["standard_correlation_report"] = str(correlation_path)
    if date_tag is not None:
        payload["standard_date_tag"] = date_tag
    if legacy_csv_path is not None:
        payload["legacy_output"] = str(legacy_csv_path)
    if extras_payload:
        payload.update(extras_payload)
    if metrics is not None:
        summary = metrics.summary()
        if summary.get("rows") is not None:
            payload["postprocess_rows"] = summary["rows"]
        if summary.get("columns") is not None:
            payload["postprocess_columns"] = summary["columns"]
        if summary.get("duration_s") is not None:
            payload["postprocess_duration_s"] = summary["duration_s"]
        if summary.get("steps") is not None:
            payload["postprocess_steps"] = summary["steps"]
        validation = getattr(metrics, "validation", None)
        if validation is not None and getattr(validation, "schema", None) is not None:
            payload["postprocess_schema"] = validation.schema
    logger.info(event, **payload)


def run_pubmed(
    cfg: Config,
    args: argparse.Namespace,
    *,
    pipeline: DocumentPipeline | None = None,
) -> int:
    """Execute the ``pubmed`` mode.

    Parameters
    ----------
    cfg : Config
        Application configuration containing rate limiting, API and CSV export
        settings.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. A non-zero value indicates that an error occurred
        while reading input identifiers, fetching metadata or writing the
        resulting CSV.
    """
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
            args.output_csv = output_path
        else:
            output_path = Path(
                io.default_output_path(
                    args.input_csv,
                    cfg.io,
                    date=getattr(args, "date", None),
                )
            )
            args.final_out = output_path
            args.output_csv = output_path
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        args.output_csv = output_path
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
            exc_info=exc,
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

    partial_run = (limit is not None) or (offset > 0)

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
        delimiter = getattr(args, "fallback_doi_delimiter", None) or cfg.io.csv_sep
        encoding = getattr(args, "fallback_doi_encoding", None) or cfg.io.csv_encoding
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
                exc_info=exc,
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
                exc_info=exc,
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
            fallback_doi_map=(fallback_state.mapping if fallback_state else None),
            return_generator=True,
        )
        output = output_path
        if fallback_state is not None:
            fallback_state.metrics.mark_total_candidates(len(fallback_state.mapping))

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
        finalise_result = _finalise_export(
            normalised_frames,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
            chunk_size=getattr(args, "batch_size", pubmed_defaults.batch_size),
            partial_run=partial_run,
            rerun_postprocess=bool(getattr(args, "rerun_postprocess", False)),
            postprocess_enabled=bool(getattr(args, "postprocess", False)),
            postprocess_config_path=getattr(args, "config", None),
            emit_legacy_artifacts=bool(getattr(args, "emit_legacy_artifacts", False)),
            date_tag=getattr(args, "_standard_date_tag", None),
            cli_args=args,
            metadata_sources=_PUBMED_METADATA_SOURCES,
            stats_extra=service.stats_extra,
        )
        exit_code = finalise_result.exit_code
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error(
            "pubmed_pipeline_failed",
            error=str(exc),
            exc_info=exc,
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
    output_dataset = (
        finalise_result.artifacts.dataset
        if finalise_result.artifacts is not None
        else output_path
    )
    if exit_code == 0:
        _log_document_completion(
            "document_pubmed_done",
            cfg=cfg,
            output_path=output_dataset,
            logger=logger,
            extras={"mode": "pubmed"},
            finalise_result=finalise_result,
            postprocess_enabled=bool(getattr(args, "postprocess", False)),
        )
    else:
        logger.error(
            "document_pubmed_failed",
            output=str(output_dataset),
            exit_code=exit_code,
        )
    return exit_code


def run_chembl(
    cfg: Config,
    args: argparse.Namespace,
    *,
    pipeline: DocumentPipeline | None = None,
) -> int:
    """Execute the ``chembl`` mode.

    Parameters
    ----------
    cfg : Config
        Application configuration providing ChEMBL client, retry and CSV export
        options.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. A non-zero value indicates that reading the input
        identifiers, fetching ChEMBL data or exporting the CSV failed. Network
        errors are logged and converted into placeholder rows where possible.
    """
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
            args.output_csv = output_path
        else:
            output_path = Path(
                io.default_output_path(
                    args.input_csv,
                    cfg.io,
                    date=getattr(args, "date", None),
                )
            )
            args.final_out = output_path
            args.output_csv = output_path
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        args.output_csv = output_path
    chunk_size = getattr(args, "chunk_size", chembl_defaults.chunk_size)
    timeout = _resolve_timeout(getattr(args, "timeout", None), chembl_defaults.timeout)
    args.timeout = timeout
    metadata_obj = getattr(args, "_config_metadata", None)
    if not isinstance(metadata_obj, ConfigMetadata):
        metadata_obj = None
    output_source = "cli" if getattr(args, "final_out", None) else "derived"
    service = pipeline or DocumentPipeline(cfg)
    logger.info(
        "document_chembl_start",
        input=prepare_option(
            metadata_obj, value=str(args.input_csv), default_source="cli"
        ),
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

    # Configure session for ChEMBL requests
    with ETLContext(cfg) as context:
        client = context.chembl_client
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
                exc_info=exc,
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

        partial_run = (limit is not None) or (offset > 0)
        try:
            df = cl.get_documents(
                ids,
                cfg=cfg.api,
                client=client,
                chunk_size=getattr(args, "chunk_size", chembl_defaults.chunk_size),
                timeout=timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error(
                "chembl_documents_fetch_failed",
                error=str(exc),
                exc_info=exc,
                chunk_size=getattr(args, "chunk_size", chembl_defaults.chunk_size),
                timeout=timeout,
            )
            return 1
        if "doi" in df.columns:
            df["doi"] = df["doi"].map(normalise_doi)
        output = output_path
        df = normalize_documents(df)
        finalise_result = _finalise_export(
            df,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
            chunk_size=getattr(args, "chunk_size", chembl_defaults.chunk_size),
            partial_run=partial_run,
            rerun_postprocess=bool(getattr(args, "rerun_postprocess", False)),
            postprocess_enabled=bool(getattr(args, "postprocess", False)),
            postprocess_config_path=getattr(args, "config", None),
            emit_legacy_artifacts=bool(getattr(args, "emit_legacy_artifacts", False)),
            date_tag=getattr(args, "_standard_date_tag", None),
            cli_args=args,
            metadata_sources=_CHEMBL_METADATA_SOURCES,
        )
        exit_code = finalise_result.exit_code
        if exit_code == 0:
            extras: dict[str, object] = {"mode": "chembl"}
            if partial_run:
                extras["partial_run"] = True
            output_dataset = (
                finalise_result.artifacts.dataset
                if finalise_result.artifacts is not None
                else output_path
            )
            _log_document_completion(
                "document_chembl_done",
                cfg=cfg,
                output_path=output_dataset,
                logger=logger,
                extras=extras,
                finalise_result=finalise_result,
                postprocess_enabled=bool(getattr(args, "postprocess", False)),
            )
        else:
            logger.error(
                "document_chembl_failed",
                output=str(
                    finalise_result.artifacts.dataset
                    if finalise_result.artifacts is not None
                    else output_path
                ),
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
    """Execute the ``all`` mode by merging ChEMBL and PubMed outputs.

    Parameters
    ----------
    cfg : Config
        Application configuration combining all document pipeline defaults.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. Non-zero results indicate that reading identifiers,
        fetching data from upstream APIs or writing derived artefacts failed.

    Raises
    ------
    ValueError
        Raised when DOI fallback information derived from the ChEMBL payload is
        internally inconsistent.
    """
    service = pipeline or DocumentPipeline(cfg)
    all_defaults = cfg.document.all
    limit = getattr(args, "limit", all_defaults.limit)
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="document.all", limit=limit)
        return 1

    # Prepare shared session before performing any API calls
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
            exc_info=exc,
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

    partial_run = (limit is not None) or (offset > 0)
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
            args.output_csv = output_path
        else:
            output_path = Path(
                io.default_output_path(
                    args.input_csv,
                    cfg.io,
                    date=getattr(args, "date", None),
                )
            )
            args.final_out = output_path
            args.output_csv = output_path
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        args.output_csv = output_path
    fallback_enabled = getattr(args, "fallback_doi_enabled", False)
    fallback_path_arg = getattr(args, "fallback_doi_path", None)
    chembl_chunk_size = getattr(args, "chembl_chunk_size", all_defaults.chunk_size)
    args.chembl_chunk_size = chembl_chunk_size
    chembl_timeout_value = getattr(args, "chembl_timeout", None)
    if chembl_timeout_value is None:
        chembl_timeout_value = getattr(args, "timeout", None)
    chembl_timeout = _resolve_timeout(chembl_timeout_value, all_defaults.timeout)
    args.chembl_timeout = chembl_timeout
    pubmed_timeout_value = getattr(args, "pubmed_timeout", None)
    if pubmed_timeout_value is None:
        pubmed_timeout_value = getattr(args, "timeout", None)
    pubmed_timeout = _resolve_timeout(pubmed_timeout_value, PUBMED_DEFAULTS.timeout)
    args.pubmed_timeout = pubmed_timeout
    logger.info(
        "document_all_start",
        input=str(args.input_csv),
        output=str(output_path),
        limit=limit,
        offset=offset,
        chembl_chunk_size=chembl_chunk_size,
        pubmed_workers=getattr(args, "pubmed_workers", all_defaults.workers),
        pubmed_batch_size=getattr(args, "pubmed_batch_size", all_defaults.batch_size),
        pubmed_sleep=getattr(args, "pubmed_sleep", all_defaults.sleep),
        chembl_timeout=chembl_timeout,
        pubmed_timeout=pubmed_timeout,
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
        delimiter = getattr(args, "fallback_doi_delimiter", None) or cfg.io.csv_sep
        encoding = getattr(args, "fallback_doi_encoding", None) or cfg.io.csv_encoding
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
                exc_info=exc,
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
                exc_info=exc,
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
        with ETLContext(cfg) as context:
            client = context.chembl_client
            doc_df = cl.get_documents(
                ids_for_fetch,
                cfg=cfg.api,
                client=client,
                chunk_size=chembl_chunk_size,
                timeout=chembl_timeout,
            )
    except (requests.RequestException, ValueError) as exc:
        logger.error(
            "chembl_documents_fetch_failed",
            ids=sample_ids,
            error=str(exc),
            exc_info=exc,
            output=str(output_path),
            chunk_size=chembl_chunk_size,
            timeout=chembl_timeout,
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
        finalise_result = _finalise_export(
            processed,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
            chunk_size=getattr(args, "chembl_chunk_size", all_defaults.chunk_size),
            partial_run=partial_run,
            rerun_postprocess=bool(getattr(args, "rerun_postprocess", False)),
            postprocess_enabled=bool(getattr(args, "postprocess", False)),
            postprocess_config_path=getattr(args, "config", None),
            emit_legacy_artifacts=bool(getattr(args, "emit_legacy_artifacts", False)),
            date_tag=getattr(args, "_standard_date_tag", None),
            cli_args=args,
            metadata_sources=_ALL_METADATA_SOURCES,
        )
        exit_code = finalise_result.exit_code
        if fallback_state is not None:
            logger.info(
                "fallback_doi_metrics",
                pipeline="all",
                path=str(fallback_state.path),
                overwrite=fallback_state.overwrite,
                **fallback_state.metrics.as_log_kwargs(),
            )
        if exit_code == 0:
            extras: dict[str, object] = {"mode": "all"}
            if partial_run:
                extras["partial_run"] = True
            output_dataset = (
                finalise_result.artifacts.dataset
                if finalise_result.artifacts is not None
                else output_path
            )
            _log_document_completion(
                "document_all_done",
                cfg=cfg,
                output_path=output_dataset,
                logger=logger,
                extras=extras,
                finalise_result=finalise_result,
                postprocess_enabled=bool(getattr(args, "postprocess", False)),
            )
        else:
            logger.error(
                "document_all_failed",
                output=str(
                    finalise_result.artifacts.dataset
                    if finalise_result.artifacts is not None
                    else output_path
                ),
                exit_code=exit_code,
            )
        return exit_code

    # Normalise PubMed identifiers to strings to avoid dtype mismatches
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
            exc_info=exc,
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
    finalise_result = _finalise_export(
        processed,
        output,
        cfg,
        input_csv=Path(args.input_csv),
        key_columns=["document_chembl_id"],
        chunk_size=getattr(args, "chembl_chunk_size", all_defaults.chunk_size),
        partial_run=partial_run,
        rerun_postprocess=bool(getattr(args, "rerun_postprocess", False)),
        postprocess_enabled=bool(getattr(args, "postprocess", False)),
        postprocess_config_path=getattr(args, "config", None),
        emit_legacy_artifacts=bool(getattr(args, "emit_legacy_artifacts", False)),
        date_tag=getattr(args, "_standard_date_tag", None),
        cli_args=args,
        metadata_sources=_ALL_METADATA_SOURCES,
        stats_extra=service.stats_extra,
    )
    exit_code = finalise_result.exit_code
    if fallback_state is not None:
        logger.info(
            "fallback_doi_metrics",
            pipeline="all",
            path=str(fallback_state.path),
            overwrite=fallback_state.overwrite,
            **fallback_state.metrics.as_log_kwargs(),
        )
    if exit_code == 0:
        extras: dict[str, object] = {"mode": "all"}
        if partial_run:
            extras["partial_run"] = True
        output_dataset = (
            finalise_result.artifacts.dataset
            if finalise_result.artifacts is not None
            else output_path
        )
        _log_document_completion(
            "document_all_done",
            cfg=cfg,
            output_path=output_dataset,
            logger=logger,
            extras=extras,
            finalise_result=finalise_result,
            postprocess_enabled=bool(getattr(args, "postprocess", False)),
        )
    else:
        logger.error(
            "document_all_failed",
            output=str(
                finalise_result.artifacts.dataset
                if finalise_result.artifacts is not None
                else output_path
            ),
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
            output_path = Path(
                io.default_output_path(
                    args.input_csv,
                    cfg.io,
                    date=getattr(args, "date", None),
                )
            )
            args.final_out = output_path
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
    args.output_csv = output_path
    table_name, inferred_date = _resolve_table_name_and_date(output_path)

    explicit_date = getattr(args, "date", None)
    if isinstance(explicit_date, str):
        explicit_date = explicit_date.strip() or None

    fallback_date = getattr(args, "date_prefix", None)
    if isinstance(fallback_date, str):
        fallback_date = fallback_date.strip() or None

    standard_date_tag = (
        explicit_date
        or fallback_date
        or inferred_date
        or datetime.now(UTC).strftime("%Y%m%d")
    )
    args._standard_date_tag = standard_date_tag
    emit_legacy = bool(getattr(args, "emit_legacy_artifacts", False))

    output_dir_value = getattr(cfg.io, "output_dir", None)
    output_dir = Path(output_dir_value) if output_dir_value else output_path.parent
    canonical_dataset = output_dir / f"output.{table_name}_{standard_date_tag}.csv"

    mode = getattr(args, "mode", None)
    if mode in (None, ""):
        mode = getattr(args, "command", None)
    chembl_timeout_override: float | None = None
    pubmed_timeout_override: float | None = None
    if mode == "chembl":
        chembl_timeout_override = getattr(args, "timeout", None)
    elif mode == "pubmed":
        pubmed_timeout_override = getattr(args, "timeout", None)
    elif mode == "all":
        chembl_timeout_override = getattr(args, "chembl_timeout", None)
        if chembl_timeout_override is None:
            chembl_timeout_override = getattr(args, "timeout", None)
        pubmed_timeout_override = getattr(args, "pubmed_timeout", None)
        if pubmed_timeout_override is None:
            pubmed_timeout_override = getattr(args, "timeout", None)
    if chembl_timeout_override is not None:
        cfg.api.timeout_read = float(chembl_timeout_override)
    if pubmed_timeout_override is not None:
        cfg.pubmed.timeout_read = float(pubmed_timeout_override)
    if args.skip_existing and not args.force:
        if emit_legacy and output_path.exists():
            logger.info("pipeline_skip_existing", output=str(output_path))
            return 0
        if (not emit_legacy) and canonical_dataset.exists():
            logger.info("pipeline_skip_existing", output=str(canonical_dataset))
            return 0
    if hasattr(args, "func") and getattr(args, "func", None) is None:
        logger.error(
            "missing_subcommand_handler",
            command=getattr(args, "command", None),
        )
        return 1
    handler = MODE_HANDLERS.get(str(mode))
    if handler is None:
        logger.error("unknown_mode", mode=mode)
        return 1
    current_handler = getattr(sys.modules[__name__], handler.__name__, handler)
    service = DocumentPipeline(cfg)

    supports_pipeline = True
    try:
        signature = inspect.signature(current_handler)
    except (TypeError, ValueError):  # pragma: no cover - builtins or C extensions
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


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the argument parser for document utilities."""

    root, _, log_cfg = build_root_parser()
    root.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
    parser = argparse.ArgumentParser(
        description="Document data utilities",
        parents=[root],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))

    if "--emit-legacy-artifacts" not in parser._option_string_actions:
        parser.add_argument(
            "--emit-legacy-artifacts",
            dest="emit_legacy_artifacts",
            action=argparse.BooleanOptionalAction,
            default=False,
            help=(
                "Write legacy CSV sidecars and metadata in addition to the standard "
                "outputs saved under io.output_dir"
            ),
        )

    set_emit_legacy_help(
        parser,
        "Write legacy CSV sidecars and metadata in addition to the standard "
        "outputs saved under io.output_dir",
    )
    parser.add_argument(
        "--rerun-postprocess",
        action="store_true",
        help=(
            "Rebuild stage-aligned exports even if a previous run already produced "
            "them"
        ),
    )
    legacy_option = parser._option_string_actions.get("--emit-legacy-artifacts")
    if legacy_option is not None:
        legacy_option.help = (
            "Write legacy CSV sidecars and metadata in addition to the standard "
            "outputs saved under io.output_dir"
        )

    pipeline_group = parser.add_argument_group("Pipeline selection")
    pipeline_group.add_argument(
        "--mode",
        choices=("chembl", "pubmed", "all"),
        required=False,
        default=None,
        help="Document pipeline to execute",
    )
    parser.add_argument(
        "command",
        nargs="?",
        choices=("chembl", "pubmed", "all"),
        help=argparse.SUPPRESS,
    )
    pipeline_group.add_argument(
        "--column",
        default=PUBMED_DEFAULTS.column,
        help=(
            "Input column containing identifiers (defaults: "
            f"pubmed={PUBMED_DEFAULTS.column}, chembl={CHEMBL_DEFAULTS.column}, "
            f"all={ALL_DEFAULTS.column})"
        ),
    )
    pipeline_group.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process (default: no limit)",
    )
    pipeline_group.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    pipeline_group.add_argument(
        "--postprocess",
        dest="postprocess",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Enable document postprocessing after the main pipeline",
    )
    pipeline_group.add_argument(
        "--openalex-rps",
        type=float,
        default=None,
        help="Requests per second limit for OpenAlex lookups",
    )
    pipeline_group.add_argument(
        "--crossref-rps",
        type=float,
        default=None,
        help="Requests per second limit for CrossRef lookups",
    )

    single_group = parser.add_argument_group("Single pipeline options")
    single_group.add_argument(
        "--batch-size",
        type=positive_int,
        default=PUBMED_DEFAULTS.batch_size,
        help=(
            "Maximum PMIDs per PubMed request when running in pubmed mode "
            f"(default: {PUBMED_DEFAULTS.batch_size})"
        ),
    )
    single_group.add_argument(
        "--sleep",
        type=float,
        default=PUBMED_DEFAULTS.sleep,
        help=(
            "Seconds to sleep between PubMed requests when running in pubmed mode "
            f"(default: {PUBMED_DEFAULTS.sleep})"
        ),
    )
    single_group.add_argument(
        "--workers",
        type=int,
        default=PUBMED_DEFAULTS.workers,
        help=(
            "Number of concurrent PubMed requests when running in pubmed mode "
            f"(default: {PUBMED_DEFAULTS.workers})"
        ),
    )
    single_group.add_argument(
        "--chunk-size",
        type=positive_int,
        default=CHEMBL_DEFAULTS.chunk_size,
        help=(
            "Maximum identifiers per ChEMBL request when running in chembl mode "
            f"(default: {CHEMBL_DEFAULTS.chunk_size})"
        ),
    )
    single_group.add_argument(
        "--timeout",
        type=float,
        default=None,
        help=(
            "HTTP read timeout in seconds (defaults: chembl/all="
            f"{CHEMBL_DEFAULTS.timeout}, pubmed={PUBMED_DEFAULTS.timeout})"
        ),
    )

    combined_group = parser.add_argument_group("Combined pipeline overrides")
    combined_group.add_argument(
        "--chembl-chunk-size",
        "--chembl-batch-size",
        dest="chembl_chunk_size",
        type=positive_int,
        default=ALL_DEFAULTS.chunk_size,
        help=(
            "Maximum identifiers per ChEMBL request when running in all mode "
            f"(default: {ALL_DEFAULTS.chunk_size})"
        ),
    )
    combined_group.add_argument(
        "--chembl-timeout",
        dest="chembl_timeout",
        type=float,
        default=None,
        help=(
            "Timeout in seconds for ChEMBL requests when running in all mode "
            f"(default: {ALL_DEFAULTS.timeout})"
        ),
    )
    combined_group.add_argument(
        "--pubmed-sleep",
        "--chembl-sleep",
        dest="pubmed_sleep",
        type=float,
        default=ALL_DEFAULTS.sleep,
        help=(
            "Seconds to sleep between PubMed requests when running in all mode "
            f"(default: {ALL_DEFAULTS.sleep})"
        ),
    )
    combined_group.add_argument(
        "--pubmed-workers",
        "--chembl-workers",
        dest="pubmed_workers",
        type=int,
        default=ALL_DEFAULTS.workers,
        help=(
            "Number of concurrent PubMed requests when running in all mode "
            f"(default: {ALL_DEFAULTS.workers})"
        ),
    )
    combined_group.add_argument(
        "--pubmed-batch-size",
        "--pubmed-chunk-size",
        dest="pubmed_batch_size",
        type=positive_int,
        default=ALL_DEFAULTS.batch_size,
        help=(
            "Maximum PMIDs per PubMed request when running in all mode "
            f"(default: {ALL_DEFAULTS.batch_size})"
        ),
    )
    combined_group.add_argument(
        "--pubmed-timeout",
        dest="pubmed_timeout",
        type=float,
        default=None,
        help=(
            "Timeout in seconds for PubMed requests when running in all mode "
            f"(default: {PUBMED_DEFAULTS.timeout})"
        ),
    )

    fallback_group = parser.add_argument_group("Fallback DOI overrides")
    fallback_group.add_argument(
        "--fallback-doi-enabled",
        action="store_true",
        help="Enable lookup of DOI overrides from a CSV file",
    )
    parser.set_defaults(fallback_doi_csv=None)
    fallback_group.add_argument(
        "--fallback-doi-path",
        "--fallback-doi-csv",
        dest="fallback_doi_path",
        action=_FallbackPathAction,
        type=path_argument,
        default=None,
        help="CSV file containing DOI overrides keyed by PMID",
    )
    fallback_group.add_argument(
        "--fallback-doi-col-pmid",
        default="PMID",
        help="Column containing PubMed identifiers in the fallback CSV",
    )
    fallback_group.add_argument(
        "--fallback-doi-col-doi",
        default="DOI",
        help="Column containing DOI values in the fallback CSV",
    )
    fallback_group.add_argument(
        "--fallback-doi-delimiter",
        default=None,
        help="Delimiter used when reading the fallback CSV (default: io.csv_sep)",
    )
    fallback_group.add_argument(
        "--fallback-doi-encoding",
        default=None,
        help="Encoding used for the fallback CSV (default: io.csv_encoding)",
    )
    fallback_group.add_argument(
        "--fallback-doi-overwrite",
        action="store_true",
        help="Allow replacing existing DOIs with fallback values",
    )

    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults.

    Parameters
    ----------
    argv : Sequence[str] | None, optional
        Command-line arguments to parse. When ``None`` the process arguments
        from :data:`sys.argv` are used.

    Returns
    -------
    int
        ``0`` when the selected pipeline completes successfully, non-zero
        otherwise.

    Raises
    ------
    SystemExit
        Raised when argument validation fails and ``argparse`` aborts
        execution.
    """
    parser, log_cfg = build_parser()
    if argv is None:
        argv_list: list[str] = list(sys.argv[1:])
    else:
        argv_list = [str(item) for item in argv]
    args = parser.parse_args(argv_list)
    prepare_io_paths(
        args,
        input_default=DEFAULT_INPUT_NAME,
        output_stem=DEFAULT_OUTPUT_STEM,
    )
    limit_value = getattr(args, "limit", None)
    if limit_value == 0:
        logger.info("pipeline_skip_limit", limit=limit_value)
        return 0
    mode = getattr(args, "mode", None)
    command_value = getattr(args, "command", None)
    if not mode and command_value:
        mode = command_value
        args.mode = mode
    if not mode:
        parser.error("--mode is required")
    if command_value is None and mode is not None:
        args.command = mode
    mode = str(mode)
    if limit_value is not None and limit_value < 0:
        parser.error("--limit must be zero or a positive integer")
    offset_value = getattr(args, "offset", 0)
    if offset_value < 0:
        parser.error("--offset must be zero or a positive integer")
    if mode in {"pubmed", "all"}:
        fallback_enabled_cli = getattr(args, "fallback_doi_enabled", False)
        fallback_path_cli = getattr(args, "fallback_doi_path", None)
        if fallback_path_cli is not None and not fallback_enabled_cli:
            fallback_enabled_cli = True
            args.fallback_doi_enabled = True
        if fallback_enabled_cli:
            if fallback_path_cli is None:
                parser.error(
                    "--fallback-doi-path is required when fallback DOI overrides are enabled"
                )
            fallback_path = (
                fallback_path_cli
                if isinstance(fallback_path_cli, Path)
                else Path(str(fallback_path_cli))
            )
            if not fallback_path.exists():
                parser.error("--fallback-doi-path must point to an existing file")
            if not fallback_path.is_file():
                parser.error("--fallback-doi-path must be a file")
            if not os.access(fallback_path, os.R_OK):
                parser.error("--fallback-doi-path must be readable")
            delimiter_cli = getattr(args, "fallback_doi_delimiter", None)
            if delimiter_cli is not None:
                delimiter_text = str(delimiter_cli)
                if not delimiter_text:
                    parser.error("--fallback-doi-delimiter must not be empty")
                if len(delimiter_text) > 1:
                    parser.error("--fallback-doi-delimiter must be a single character")
            encoding_cli = getattr(args, "fallback_doi_encoding", None)
            if encoding_cli is not None and not str(encoding_cli).strip():
                parser.error("--fallback-doi-encoding must not be empty")
            for attr_name in ("fallback_doi_col_pmid", "fallback_doi_col_doi"):
                attr_value = getattr(args, attr_name, None)
                if attr_value is None or not str(attr_value).strip():
                    option = attr_name.replace("_", "-")
                    parser.error(f"--{option} must not be empty")
    mapping = {
        "column": f"document.{mode}.column",
        "limit": f"document.{mode}.limit",
    }
    if mode == "pubmed":
        mapping.update(
            {
                "sleep": "document.pubmed.sleep",
                "workers": "document.pubmed.workers",
                "batch_size": "document.pubmed.batch_size",
                "timeout": "sources.pubmed.timeout_read",
            }
        )
    elif mode == "chembl":
        mapping.update(
            {
                "chunk_size": "document.chembl.chunk_size",
                "timeout": "document.chembl.timeout",
            }
        )
    elif mode == "all":
        mapping.update(
            {
                "chunk_size": "document.all.chunk_size",
                "timeout": "document.all.timeout",
                "chembl_chunk_size": "document.all.chunk_size",
                "pubmed_sleep": "document.all.sleep",
                "pubmed_workers": "document.all.workers",
                "pubmed_batch_size": "document.all.batch_size",
                "chembl_timeout": "document.all.timeout",
                "pubmed_timeout": "sources.pubmed.timeout_read",
            }
        )
    mapping |= {
        "openalex_rps": "openalex.rps",
        "crossref_rps": "crossref.rps",
    }
    with setup_cli_logging(
        Path(__file__).with_suffix("").name, log_cfg, getattr(args, "date", None)
    ) as logging_ctx:
        exit_code = run_cli_command(
            args=args,
            parser=parser,
            log_cfg=logging_ctx.log_cfg,
            mapping=mapping,
            run=run,
            logger=logger,
        )
    configure_logger(log_cfg)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
