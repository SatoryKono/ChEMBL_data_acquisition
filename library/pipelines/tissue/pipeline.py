"""Orchestration helpers for the tissue data acquisition pipeline."""

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
from ...schemas import normalize_tissues
from ...validation import validate_tissues
from .chembl import TISSUE_BASE_COLUMNS, TISSUE_COLUMN_ORDER, get_tissues


@dataclass(slots=True)
class TissuePipelineOptions:
    """Configuration for executing the tissue pipeline."""

    input_csv: Path
    output_csv: Path
    column: str
    batch_size: int
    limit: int | None = None
    offset: int = 0
    timeout: float | None = None
    mode: str = "chembl"


@dataclass(slots=True)
class TissuePipelineResult:
    """Outcome of running the tissue pipeline."""

    exit_code: int
    records: int
    duration: float
    output_path: Path
    failure_path: Path | None
    failures: int
    missing_ids: tuple[str, ...]
    written: bool


def read_tissue_identifiers(
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


def _ensure_string_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """Coerce the specified ``columns`` to :class:`pandas.StringDtype`."""

    result = df.copy()
    for column in columns:
        if column not in result.columns:
            result[column] = pd.Series(pd.NA, dtype=pd.StringDtype())
        else:
            result[column] = result[column].astype(pd.StringDtype())
    return result


def prepare_tissue_dataframe(
    df: pd.DataFrame,
    requested: Sequence[str],
) -> tuple[pd.DataFrame, tuple[str, ...]]:
    """Return normalised tissue records together with missing identifiers."""

    work = df.copy()
    for column in TISSUE_BASE_COLUMNS:
        if column not in work.columns:
            work[column] = pd.NA
    work = work.reindex(columns=TISSUE_BASE_COLUMNS, fill_value=pd.NA)
    work = _ensure_string_columns(work, TISSUE_BASE_COLUMNS)
    work = work.replace({"": pd.NA})

    present = {
        value.strip().upper()
        for value in work["tissue_chembl_id"].dropna()
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
                for column in TISSUE_BASE_COLUMNS
            }
        )
        placeholder["tissue_chembl_id"] = pd.Series(missing, dtype=pd.StringDtype())
        work = pd.concat([work, placeholder], ignore_index=True, sort=False)

    work = work.drop_duplicates(subset=["tissue_chembl_id"], keep="first")
    work = work.reindex(columns=TISSUE_BASE_COLUMNS, fill_value=pd.NA)
    work = normalize_tissues(work)
    work = _ensure_string_columns(work, TISSUE_BASE_COLUMNS)
    return work, tuple(missing)


def _write_failures(
    failure_cases: pd.DataFrame,
    *,
    destination: Path,
    cfg: Config,
) -> None:
    """Persist validation failures to ``destination`` using deterministic CSV."""

    if failure_cases.empty:
        if destination.exists():
            destination.unlink()
        return

    failure_key_cols = [
        col for col in ("index", "column") if col in failure_cases.columns
    ]
    if not failure_key_cols and not failure_cases.empty:
        failure_key_cols = [failure_cases.columns[0]]

    write_csv_deterministic(
        failure_cases.copy(),
        destination,
        col_order=list(failure_cases.columns),
        key_cols=failure_key_cols,
        sep=cfg.io.csv_sep,
        encoding=cfg.io.csv_encoding,
        chunksize=cfg.io.csv_chunksize,
        sort_chunksize=cfg.io.csv_chunksize,
        cfg=cfg,
    )


def run_tissue_pipeline(
    cfg: Config,
    options: TissuePipelineOptions,
    *,
    client: ChemblClient | None = None,
) -> TissuePipelineResult:
    """Execute the tissue retrieval pipeline and return summary information."""

    start = perf_counter()
    failure_path = options.output_csv.with_name(
        f"{options.output_csv.stem}_validation_failures.csv"
    )
    identifiers: list[str]
    try:
        identifiers = read_tissue_identifiers(
            options.input_csv,
            column=options.column,
            io_cfg=cfg.io,
            limit=options.limit,
            offset=options.offset,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("read_fail", error=str(exc), path=str(options.input_csv))
        duration = perf_counter() - start
        return TissuePipelineResult(
            exit_code=1,
            records=0,
            duration=duration,
            output_path=options.output_csv,
            failure_path=None,
            failures=0,
            missing_ids=(),
            written=False,
        )

    logger.info("tissue_identifiers_loaded", count=len(identifiers))

    if not identifiers:
        empty = pd.DataFrame(columns=TISSUE_BASE_COLUMNS)
        empty = add_pipeline_metadata(empty)
        empty = normalize_tissues(empty)
        empty = _ensure_string_columns(empty, TISSUE_COLUMN_ORDER)
        write_csv_deterministic(
            empty.copy(),
            options.output_csv,
            col_order=TISSUE_COLUMN_ORDER,
            key_cols=["tissue_chembl_id"],
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            chunksize=cfg.io.csv_chunksize,
            sort_chunksize=cfg.io.csv_chunksize,
            cfg=cfg,
        )
        duration = perf_counter() - start
        return TissuePipelineResult(
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
        logger.info("tissue_fetch_start", requested=len(identifiers))
        effective_timeout = (
            options.timeout if options.timeout is not None else cfg.tissue.timeout
        )
        raw = get_tissues(
            identifiers,
            cfg=cfg.api,
            client=session_client,
            chunk_size=max(1, options.batch_size),
            timeout=effective_timeout,
        )
    except (requests.RequestException, ValueError) as exc:
        logger.error("tissue_fetch_failed", error=str(exc))
        duration = perf_counter() - start
        if close_client:
            session_client.close()
        return TissuePipelineResult(
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

    prepared, missing_ids = prepare_tissue_dataframe(raw, identifiers)
    if missing_ids:
        sample = list(missing_ids[:5])
        logger.warning(
            "tissue_missing_identifiers",
            total=len(missing_ids),
            sample=sample,
        )

    enriched = add_pipeline_metadata(prepared)
    enriched = _ensure_string_columns(enriched, TISSUE_COLUMN_ORDER)
    enriched = normalize_tissues(enriched)
    enriched = _ensure_string_columns(enriched, TISSUE_COLUMN_ORDER)
    enriched = enriched.reindex(columns=TISSUE_COLUMN_ORDER, fill_value=pd.NA)

    validation = validate_tissues(enriched, return_result=True)
    failure_cases = validation.failure_cases
    if not failure_cases.empty:
        _write_failures(failure_cases, destination=failure_path, cfg=cfg)
        duration = perf_counter() - start
        return TissuePipelineResult(
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
    output_df = _ensure_string_columns(output_df, TISSUE_COLUMN_ORDER)
    output_df = output_df.reindex(columns=TISSUE_COLUMN_ORDER, fill_value=pd.NA)

    write_csv_deterministic(
        output_df,
        options.output_csv,
        col_order=TISSUE_COLUMN_ORDER,
        key_cols=["tissue_chembl_id"],
        sep=cfg.io.csv_sep,
        encoding=cfg.io.csv_encoding,
        chunksize=cfg.io.csv_chunksize,
        sort_chunksize=cfg.io.csv_chunksize,
        cfg=cfg,
    )

    duration = perf_counter() - start
    return TissuePipelineResult(
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
    "TissuePipelineOptions",
    "TissuePipelineResult",
    "prepare_tissue_dataframe",
    "read_tissue_identifiers",
    "run_tissue_pipeline",
]
