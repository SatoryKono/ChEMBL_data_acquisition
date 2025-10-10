"""Orchestration helpers for the cell line data acquisition pipeline."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from time import perf_counter

import pandas as pd
import requests

from ... import io
from ...clients import ChemblClient
from ...common.csv_utils import write_csv_deterministic
from ...common.log import logger
from ...common.rate_limiter import get_global_limiter
from ...config import Config, IoCfg
from ...pipelines.common import add_pipeline_metadata
from ...schemas import normalize_cell_lines
from ...validation import validate_celllines
from .chembl import (
    CELL_LINE_BASE_COLUMNS,
    CELL_LINE_COLUMN_ORDER,
    get_cell_lines,
)


@dataclass(slots=True)
class CellLinePipelineOptions:
    """Configuration for executing the cell line pipeline."""

    input_csv: Path
    output_csv: Path
    column: str
    batch_size: int
    limit: int | None = None
    offset: int = 0
    timeout: float | None = None
    mode: str = "chembl"


@dataclass(slots=True)
class CellLinePipelineResult:
    """Outcome of running the cell line pipeline."""

    exit_code: int
    records: int
    duration: float
    output_path: Path
    failure_path: Path | None
    failures: int
    missing_ids: tuple[str, ...]
    written: bool


def read_cellline_identifiers(
    input_csv: Path,
    *,
    column: str,
    io_cfg: IoCfg,
    limit: int | None,
    offset: int,
) -> list[str]:
    """Return identifiers from ``input_csv`` honouring ``limit`` and ``offset``."""

    ids_iter = io.read_ids(input_csv, column=column, cfg=io_cfg)
    ids_iter = (identifier.strip() for identifier in ids_iter)
    if offset:
        ids_iter = islice(ids_iter, offset, None)
    if limit is not None:
        ids_iter = islice(ids_iter, limit)
    identifiers = [identifier for identifier in ids_iter if identifier]
    return identifiers


def _ensure_column_types(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce cell line columns to their canonical nullable dtypes."""

    result = df.copy()
    string_columns = [
        "cell_chembl_id",
        "cell_name",
        "cell_description",
        "cell_source_organism",
        "cell_source_tissue",
        "cellosaurus_id",
        "cl_lincs_id",
        "clo_id",
        "efo_id",
    ]
    integer_columns = ["cell_id", "cell_source_tax_id"]

    for column in CELL_LINE_BASE_COLUMNS:
        if column not in result.columns:
            result[column] = pd.NA

    for column in string_columns:
        result[column] = result[column].astype(pd.StringDtype())

    for column in integer_columns:
        numeric = pd.to_numeric(result[column], errors="coerce")
        result[column] = numeric.astype("Int64")

    return result


def prepare_cellline_dataframe(
    df: pd.DataFrame,
    requested: Sequence[str],
) -> tuple[pd.DataFrame, tuple[str, ...]]:
    """Return normalised cell line records together with missing identifiers."""

    work = df.copy()
    for column in CELL_LINE_BASE_COLUMNS:
        if column not in work.columns:
            work[column] = pd.NA
    work = work.reindex(columns=CELL_LINE_BASE_COLUMNS, fill_value=pd.NA)
    work = work.replace({"": pd.NA})
    work = _ensure_column_types(work)

    present = {
        value.strip().upper()
        for value in work["cell_chembl_id"].dropna()
        if isinstance(value, str) and value.strip()
    }
    missing: list[str] = []
    for identifier in requested:
        normalised = str(identifier).strip()
        if not normalised:
            continue
        if normalised.upper() not in present:
            missing.append(identifier)

    if missing:
        placeholder = pd.DataFrame(
            {
                column: pd.Series(
                    pd.NA, index=range(len(missing)), dtype=pd.StringDtype()
                )
                for column in CELL_LINE_BASE_COLUMNS
            }
        )
        placeholder["cell_chembl_id"] = pd.Series(missing, dtype=pd.StringDtype())
        work = pd.concat([work, placeholder], ignore_index=True, sort=False)

    work = work.drop_duplicates(subset=["cell_chembl_id"], keep="first")
    work = work.reindex(columns=CELL_LINE_BASE_COLUMNS, fill_value=pd.NA)
    work = _ensure_column_types(work)
    return work, tuple(missing)


def run_cellline_pipeline(
    cfg: Config,
    options: CellLinePipelineOptions,
    *,
    client: ChemblClient | None = None,
) -> CellLinePipelineResult:
    """Execute the cell line retrieval pipeline and return summary information."""

    start = perf_counter()
    failure_path = options.output_csv.with_name(
        f"{options.output_csv.stem}_validation_failures.csv"
    )
    identifiers: list[str]
    try:
        identifiers = read_cellline_identifiers(
            options.input_csv,
            column=options.column,
            io_cfg=cfg.io,
            limit=options.limit,
            offset=options.offset,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("cellline_read_fail", error=str(exc), path=str(options.input_csv))
        duration = perf_counter() - start
        return CellLinePipelineResult(
            exit_code=1,
            records=0,
            duration=duration,
            output_path=options.output_csv,
            failure_path=None,
            failures=0,
            missing_ids=(),
            written=False,
        )

    logger.info("cellline_identifiers_loaded", count=len(identifiers))

    if not identifiers:
        empty = pd.DataFrame(columns=CELL_LINE_BASE_COLUMNS)
        empty = add_pipeline_metadata(empty)
        empty = normalize_cell_lines(empty)
        empty = empty.reindex(columns=CELL_LINE_COLUMN_ORDER, fill_value=pd.NA)
        empty = _ensure_column_types(empty)
        write_csv_deterministic(
            empty.copy(),
            options.output_csv,
            col_order=CELL_LINE_COLUMN_ORDER,
            key_cols=["cell_chembl_id"],
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            chunksize=cfg.io.csv_chunksize,
            sort_chunksize=cfg.io.csv_chunksize,
            cfg=cfg,
        )
        duration = perf_counter() - start
        return CellLinePipelineResult(
            exit_code=0,
            records=0,
            duration=duration,
            output_path=options.output_csv,
            failure_path=None,
            failures=0,
            missing_ids=(),
            written=True,
        )

    rate_cfg = cfg.rate
    global_limiter = None
    if (rate_cfg.global_rps or 0) > 0:
        global_limiter = get_global_limiter(rate_cfg.global_rps, rate_cfg.global_burst)

    session_client = client
    close_client = False
    if session_client is None:
        session_client = ChemblClient(
            cfg.api,
            cfg.retry,
            cfg.chembl,
            global_limiter=global_limiter,
        )
        close_client = True

    try:
        logger.info("cellline_fetch_start", requested=len(identifiers))
        effective_timeout = (
            options.timeout if options.timeout is not None else cfg.cellline.timeout
        )
        raw = get_cell_lines(
            identifiers,
            cfg=cfg.api,
            client=session_client,
            chunk_size=max(1, options.batch_size),
            timeout=effective_timeout,
        )
    except (requests.RequestException, ValueError) as exc:
        logger.error("cellline_fetch_failed", error=str(exc))
        duration = perf_counter() - start
        if close_client:
            session_client.close()
        return CellLinePipelineResult(
            exit_code=1,
            records=0,
            duration=duration,
            output_path=options.output_csv,
            failure_path=None,
            failures=0,
            missing_ids=(),
            written=False,
        )
    finally:
        if close_client:
            session_client.close()

    prepared, missing_ids = prepare_cellline_dataframe(raw, identifiers)
    if missing_ids:
        sample = list(missing_ids[:5])
        logger.warning(
            "cellline_missing_identifiers",
            total=len(missing_ids),
            sample=sample,
        )

    enriched = add_pipeline_metadata(prepared)
    enriched = _ensure_column_types(enriched)
    enriched = normalize_cell_lines(enriched)
    enriched = _ensure_column_types(enriched)
    enriched = enriched.reindex(columns=CELL_LINE_COLUMN_ORDER, fill_value=pd.NA)

    validation = validate_celllines(enriched, return_result=True)
    failure_cases = validation.failure_cases
    if not failure_cases.empty:
        write_csv_deterministic(
            failure_cases.copy(),
            failure_path,
            col_order=list(failure_cases.columns),
            key_cols=[
                col for col in ("index", "column") if col in failure_cases.columns
            ],
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            chunksize=cfg.io.csv_chunksize,
            sort_chunksize=cfg.io.csv_chunksize,
            cfg=cfg,
        )
        duration = perf_counter() - start
        return CellLinePipelineResult(
            exit_code=1,
            records=len(validation.data),
            duration=duration,
            output_path=options.output_csv,
            failure_path=failure_path,
            failures=len(failure_cases),
            missing_ids=missing_ids,
            written=False,
        )

    if failure_path.exists():
        failure_path.unlink()

    output_df = validation.data.copy()
    output_df = _ensure_column_types(output_df)
    output_df = output_df.reindex(columns=CELL_LINE_COLUMN_ORDER, fill_value=pd.NA)

    write_csv_deterministic(
        output_df,
        options.output_csv,
        col_order=CELL_LINE_COLUMN_ORDER,
        key_cols=["cell_chembl_id"],
        sep=cfg.io.csv_sep,
        encoding=cfg.io.csv_encoding,
        chunksize=cfg.io.csv_chunksize,
        sort_chunksize=cfg.io.csv_chunksize,
        cfg=cfg,
    )

    duration = perf_counter() - start
    return CellLinePipelineResult(
        exit_code=0,
        records=len(output_df),
        duration=duration,
        output_path=options.output_csv,
        failure_path=None,
        failures=0,
        missing_ids=missing_ids,
        written=True,
    )


__all__ = [
    "CellLinePipelineOptions",
    "CellLinePipelineResult",
    "prepare_cellline_dataframe",
    "read_cellline_identifiers",
    "run_cellline_pipeline",
]
